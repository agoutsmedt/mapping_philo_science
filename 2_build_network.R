#' ---
#' title: "Script for building the different types of networks"
#' author: "Aurélien Goutsmedt"
#' date: "/ Last compiled on `r format(Sys.Date())`"
#' output: 
#'   github_document:
#'     toc: true
#'     number_sections: true
#' ---

#+ r setup, include = FALSE
knitr::opts_chunk$set(eval = FALSE,
                      message = FALSE,
                      warning = FALSE)

#' # What is this script for?
#' 
#' In this script, we build different networks (cocitation, coupling and semantic similarity)
#' for different subperiods. We use the Leiden algorithm to identify communities,
#' and calculate spatial coordinates using Force Atlas 2 algorithm.
#' 
#' > WARNING: This script represents a first step of the project, and some processes have been
#' improved (notably by the creation of new functions).
#' 
#' # Loading packages, paths and data
#' 
#' ## External scripts

#+ r load_scripts, eval = TRUE
source("functions.R")
source("packages_and_paths.R")

#' ## Loading Data

#+ r load_files, eval = TRUE
direct_citations <- readRDS(paste0(data_path, "cleaned_direct_citations.rds")) %>% 
  as.data.table()
citing_articles <- readRDS(paste0(data_path, "cleaned_citing_articles.rds")) %>% 
  as.data.table()

#' # Bibliographic Coupling network
#' 
#' The first step is to build the network from the citations data. Then we will 
#' identify different communities through the Leiden algorithm, and use the 
#' Force Atlas 2 algorithm to find spatial coordinates for each citing document.
#' Finally we will be able to project the graph.
#' 
#' ## Building the network
#' 

#+ r threshold, eval = TRUE
edges_threshold <- 2
nodes_coupling_threshold <- 0

#' We have to make first choices here. We have to decide if we take all the citing
#' articles, even if they are not cited at all by other documents in our corpus.
#' Here, we have decided to 
#' `r if(nodes_coupling_threshold == 0){paste(" include all citing documents.")} else {paste(" include documents that are cited at least ", nodes_coupling_threshold, " time by other documents of the corpus.")}`
#' For the links between articles, we use a simple threshold method: we consider
#' a link between two articles important only if they share `r edges_threshold` references
#' in common.
#' 
#' We have all the authors listed in the `citing_articles` table, so we keep only the first authors for the
#' coupling. We also keep only the needed information.
#' 
#' 

#+ r keep_first_author, eval = TRUE

nodes_coupling <- citing_articles[Seq_No == 1, c("OST_BK","new_ID_Ref","Annee_Bibliographique","Full_Name","revue","titre","Abstract")]

#+ r calculate_citations, eval = TRUE
nb_cit <- direct_citations[, nb_cit := .N, by = "new_ID_Ref"]
nodes_coupling <- merge(nodes_coupling, unique(nb_cit[, c("new_ID_Ref", "nb_cit")]), by = "new_ID_Ref", all.x = "TRUE")

rm(nb_cit) # not useful anymore
# replacing NA value by 0 in nb_cit (it means that the nodes will be used in Leiden and Force Atlas, but
# not displayed in the graph)
nodes_coupling[is.na(nb_cit), ]$nb_cit <- 0

#' We will look into the distribution of citations in our corpus and will
#' reduce the number of nodes depending on the threshold (`nodes_coupling_threshold`)
#' we have chosen.
#' 

#+ r distrib_citations, eval = TRUE, message = TRUE
ggplot(nodes_coupling, aes(x = nb_cit)) +
  geom_boxplot(fill = "blue")

# reducing the number of nodes depending of the number of citations
nodes_coupling <- nodes_coupling[nb_cit >= nodes_coupling_threshold]

#' The next step just aim at securing that we only take citing documents in the
#' `direct_citations` that are also within the nodes. We also remove articles
#' without any references

#+ r build_network, eval = TRUE
edges_coupling <- unique(direct_citations[OST_BK %in% nodes_coupling$OST_BK])
if(length(edges_coupling$OST_BK) == length(direct_citations$OST_BK)){
  message("Everything was ok in the data extraction from the OST db: all the citing documents in the 
          direct citations table are also in the citing documents table.")
}
nodes_coupling <- nodes_coupling[OST_BK %in% edges_coupling$OST_BK]

#' We can now create the list of edges. We will test three different methods (see the
#' [biblionetwork](https://agoutsmedt.github.io/biblionetwork/) package documentation for more information):
#' 
#' - the coupling angle value
#' - the coupling strength value 
#' - the coupling similarity value.
#' 
#' We will check the distribution of weights for each method and choose what could be the most appropriate one.

#+ r test_edges_distribution, eval = TRUE

plot_distrib <- weight_distrib(edges_coupling, rescale = TRUE, source = "OST_BK", ref = "new_ID_Ref", weight_threshold = edges_threshold)
plot_distrib

#' Regarding the distributions of the coupling strength and the coupling similarity,
#' which are much less widespread, it seems preferable, in a first time to opt for 
#' the coupling angle measure, which will allow more difference in the weights of edges.

#+ r build_network_3, eval = TRUE
edges_coupling_angle <- biblio_coupling(edges_coupling, source = "OST_BK", ref = "new_ID_Ref", weight_threshold = edges_threshold)

#' We remove the nodes that are not any more linked because they did not pass the threshold
#' with any other articles. We can now use the tidygraph and the networkflow packages
#' and build the two networks.

#+ r build_network_4, eval = TRUE
nodes_coupling_angle <- nodes_coupling[OST_BK %in% c(edges_coupling_angle$from,edges_coupling_angle$to)]
nodes_coupling_angle[, new_ID_Ref := as.character(OST_BK)] # necessary step for tidygraph

graph_coupling <- tbl_main_component(edges_coupling_angle, nodes_coupling_angle, directed = FALSE, node_key = "OST_BK")

#' With the `edges_threshold` fixed at `r edges_threshold`, we have removed
#' `r length(nodes_coupling$OST_BK) - length(nodes_coupling_angle$OST_BK)` documents.

#' ## Clustering
#' 
#' We will now identify different communities in the network using 
#' [Leiden](https://www.nature.com/articles/s41598-019-41695-z) algorithm. 
#' 

#+ r plot_partition, eval = TRUE

plot_partition <- compare_partition(graph = graph_coupling)
plot_partition

#' Here, a resolution of 0.8 maximises modularity while allowing a reasonable number of communities.

#+ r Leiden, eval = TRUE
set.seed(1989)
graph_coupling <- leiden_workflow(graph_coupling, res_1 = 0.8)

#' We will now attribute colors to the different communities and give them a name. We give
#' two names to the communities:
#' 
#' - the ref of the most cited articles in the community
#' - the words used in the titles and abstracts of the community, with the highest tf-idf value.

#+ r naming_communities, eval = TRUE
graph_coupling <- community_colors(graph_coupling,
                                   palette = scico(n = as.integer(max(V(graph_coupling)$Com_ID)), palette = "roma"))
graph_coupling <- graph_coupling %>% 
  activate(nodes) %>% 
  mutate(words = paste(titre, Abstract),
         Label = paste0(Full_Name, "-", Annee_Bibliographique))

remove_words <- c("account",
                  "argue")

graph_coupling <- name_community_with_words(graph_coupling, title_column = "words", remove_words = remove_words)
graph_coupling <- community_names(graph_coupling, "nb_cit")

# merging the two names
graph_coupling <- graph_coupling %>% 
  mutate(Long_community_name = paste0(Community_name, "\n", com_name_by_words))

#' ## Spatialisation of the network

#+ r Layout, eval = FALSE

graph_before <- graph_coupling %>% activate(nodes) %>% mutate(id=as.character(Id),
                                                              size = nb_cit)
#nodes_before <- tbl %>% activate(nodes) %>% as.data.table()
graph_before <- graph_before %>% activate(nodes) %>% select(c(id,OST_BK,size))


write.graph(graph = graph_before, file = 'tidy.graphml', format = 'graphml')
system("java -jar GephiLayouts-1.0.jar forceatlas2 -i ./tidy.graphml -o  ./forceatlas2.graphml -threads 8 -maxiters 30000 -barneshut true -adjustsizes true -gravity 1")
gml <- read.graph("forceatlas2.graphml", format = "graphml")
graph_after <- as_tbl_graph(gml)
graph_after <- graph_after %>% activate(nodes) %>% as.data.table() %>% .[, .(x, y, OST_BK)]
graph_coupling <- graph_coupling %>% activate(nodes) %>% left_join(graph_after)

# save and load graphs
saveRDS(graph_coupling, paste0(data_path, "coupling_graph_", nodes_coupling_threshold, "-", edges_threshold,".rds"))
graph_coupling <- readRDS(paste0(data_path, "coupling_graph_", nodes_coupling_threshold, "-", edges_threshold,".rds"))

# projecting graph
top_nodes  <- top_nodes(graph_coupling, ordering_column = "nb_cit", top_n = 10, top_n_per_com = 1)
community_labels <- community_labels(graph_coupling %>% select(-Community_name) %>% rename(Community_name = Long_community_name), 
                                     community_name_column = "Community_name", community_size_column = "Size_com")

filter_nodes <- 0

agg_png(paste0(picture_path, "coupling_graph_plot_", filter_nodes, "-", nodes_coupling_threshold, "-", edges_threshold, ".png"),
        width = 60, height = 40, units = "cm", res = 300)
graph_coupling %>% 
  filter(nb_cit >= filter_nodes) %>% 
  ggraph("manual", x = x, y = y) + 
  geom_edge_arc0(aes(color = color_edges, width = weight), alpha = 0.4, strength = 0.2, show.legend = FALSE) +
  scale_edge_width_continuous(range = c(0.01,1)) +
  scale_edge_colour_identity() +
  geom_node_point(aes(x=x, y=y, size = nb_cit, fill = color), pch = 21, alpha = 0.9, show.legend = FALSE) +
  scale_size_continuous(range = c(0.01,6)) +
  scale_fill_identity() +
  new_scale("size") +
  geom_text_repel(data=top_nodes, aes(x=x, y=y, label = Label), size = 3, fontface="bold", alpha = 1, point.padding=NA, show.legend = FALSE) +
  geom_label_repel(data=community_labels, aes(x=x, y=y, label = Community_name, fill = color, size = Size_com), fontface="bold", alpha = 0.8, point.padding=NA, show.legend = FALSE) +
  scale_size_continuous(range = c(1,5)) +
  theme_void()
invisible(dev.off())
  
#' We can also try an interactive projection of the graph with sigmajs. The first
#' step is to extract the nodes and the edges, depending on the filter we have
#' chosen. Then, we have to transform nodes and edges for them to be conform
#' to the format needed by the sigmajs package.
#' 
#+ r sigmajs, eval = TRUE

graph_coupling_filtered <- graph_coupling %>% 
  filter(nb_cit >= filter_nodes) %>% 
  activate(edges) %>% 
  filter(weight >= quantile(E(graph_coupling)$weight, 0.9))

nodes <- graph_coupling_filtered %>%
  activate(nodes) %>%
  as_tibble() %>%
  rename(id = OST_BK,
         size = nb_cit) %>%
  mutate(label = paste0(Label, ", ", titre),
         id = as.character(id)) %>% 
  select(id, label, size, color, x, y) 

edges <- graph_coupling_filtered %>%
  activate(edges) %>%
  as_tibble() %>%
  rename(source = Source, 
         target = Target, 
         size = weight, 
         color = color_edges) %>%
  mutate(id = 1:n()) %>% 
  mutate(id = as.character(id),
         source = as.character(source),
         target = as.character(target)) %>% 
select(id, source, target, color, size)

# projecting the interactive graph
graph_coupling_widget <- sigmajs() %>% # initialise
  sg_nodes(nodes, id, label, size, color, x, y) %>% # add nodes
  sg_edges(edges, id, source, target, size, color) %>% # add edges
  sg_settings(drawLabels = FALSE, drawEdgeLabels = FALSE, 
              defaultEdgeType = "curve", minNodeSize = 0.1, maxNodeSize = 7, 
              minEdgeSize = 0.1, maxEdgeSize = 0.2,
              borderSize = 1, labelHoverBGColor = "node", singleHover = TRUE,
              hideEdgesOnMove = TRUE, zoomMin = 0.1, zoomMax = 1.2) %>% 
  sg_neighbors()

htmlwidgets::saveWidget(graph_coupling_widget,
                        paste0(picture_path, "coupling_graph_widget_", filter_nodes, "-", nodes_coupling_threshold, "-", edges_threshold, ".html"))

#' ## Analysing communities
#' 

tf_idf <- tf_idf(graph_coupling,
                 title_column = "words",
                 com_column = "Com_ID",
                 color_column = "color",
                 com_name_column = "Community_name",
                 com_size_column = "Size_com",
                 threshold_com = 0.01,
                 number_of_words = 12,
                 size_title_wrap = 8,
                 lemmatize_bigrams = FALSE,
                 plot = FALSE)

# plotting the graph
agg_png(paste0(picture_path, "coupling_graph_tfidf_", nodes_coupling_threshold, "-", edges_threshold, ".png"),
        width = 30, height = 20, units = "cm", res = 300)
ggplot(tf_idf[Size_com >= 0.01], aes(reorder_within(word, tf_idf, color), tf_idf, fill = color)) +
  geom_bar(stat = "identity", alpha = .8, show.legend = FALSE) +
  labs(
    title = "Highest tf-idf",
    x = "Words", y = "tf-idf"
  ) +
  facet_wrap(~Com_wrap, ncol = 4, scales = "free") +
  scale_x_reordered() +
  scale_fill_identity() +
  theme(strip.text = element_text(size = 8)) +
  coord_flip()
invisible(dev.off())
