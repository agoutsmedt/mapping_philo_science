#' ---
#' title: "List of the functions for network analysis"
#' author: "Aurélien Goutsmedt"
#' date: "/ Last compiled on `r format(Sys.Date())`"
#' output: 
#'   github_document:
#'     toc: true
#'     number_sections: true
#' ---

#+ r setup, include = FALSE
knitr::opts_chunk$set(eval = FALSE)

#' # What is this script for?
#' 
#' This script lists all the functions built for network analysis for the project. The functions are documenting in a way to be
#' later implemented in packages or at least for generating a help page (see below).
#' 
#' We first load the [docstring](https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html) package
#' to be able to write documentation in a roxygen2 style. You can then run docstring(name_of_function) to get the standard help page.

if ("docstring" %in% installed.packages() == FALSE) {
  install.packages("docstring", dependencies = TRUE)
}
library(docstring)


#' # GENERAL FUNCTIONS FOR NETWORK ANALYSIS
#' 

components_distribution <- function(edges, nodes, directed = FALSE, node_key = NULL) {
  #' Distribution of nodes in components
  #'
  #' A function which creates the graph from nodes and edges and gives a table with the percentage of nodes in each component.
  #'
  #' @param edges
  #' A dataframe with a list of links between nodes, under the columns "from" and "to". The two columns should
  #' be characters.
  #'
  #' @param nodes
  #' A dataframe with a list of nodes. The first column will be used as the identifying column.
  #' Be careful to avoid doublons in the first column. The first column should be characters.
  #'
  #' @param directed
  #' By default, the graph is computed as a non-directed graph.
  #'
  #' @param node_key
  #' The name of the identifying column in the nodes dataframe.
  #'
  #' @section Future improvements
  #' Integrate the function directly in the `tbl_main_components` function, as a second object. Imply
  #' to find how to return two objects with a function.
  
  # creating the tidygraph object
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE, node_key = node_key)
  
  # attributing a number to the different components (1 is the biggest components)
  Components_data <- graph %>%
    activate(nodes) %>%
    mutate(components_att = group_components(type = "weak")) %>%
    rename_at(1, ~"Id") %>%
    as.data.table()
  
  Components_data <- Components_data[, N_nodes_components := .N, by = "components_att"][, N_nodes_total := .N][, percentage_nodes_components := N_nodes_components / N_nodes_total * 100]
  Components_data <- unique(Components_data[order(components_att), c("components_att", "N_nodes_components", "percentage_nodes_components")])
}

weight_distrib <- function(citation_data, rescale = FALSE, ...) {
  #' Compare Weight Distribution for Different Edges Creation Methods
  #'
  #' This function creates three edges list from a data frame with citations data
  #' (a list of citing documents and cited documents). The methods used are 
  #' the [coupling angle](https://agoutsmedt.github.io/biblionetwork/reference/biblio_coupling.html),
  #' [coupling strength](https://agoutsmedt.github.io/biblionetwork/reference/coupling_strength.html)
  #' and [coupling similarity](https://agoutsmedt.github.io/biblionetwork/reference/coupling_similarity.html)
  #' measures.
  #'
  #' @param citation_data
  #' A dataframe with a list of links between nodes, under the columns "from" and "to". The two columns should
  #' be characters.
  #' 
  #' @param rescale
  #' If `TRUE`, the weights will be rescaled between 0 and 1 for each method. It
  #' makes it easier to compare the distributions of each method.
  #'
  #' @param ...
  #' Use here the parameters of the `biblio_coupling()`, `coupling_strength()`
  #' and `coupling_similarity` functions. 
 
  edges_angle <- biblionetwork::biblio_coupling(citation_data, ...)
  edges_strength <- biblionetwork::coupling_strength(citation_data, ...)
  edges_similarity <- biblionetwork::coupling_similarity(citation_data, ...)
  
  if(rescale == FALSE){
  compare_distrib <- data.table("Coupling angle" = edges_angle$weight,
                                "Coupling strength" = edges_strength$weight,
                                "Coupling similarity" = edges_similarity$weight)
  compare_distrib <- pivot_longer(compare_distrib, cols = 1:3, names_to = "Method", values_to = "Weight") %>% 
    group_by(Method) %>% 
    mutate(label = paste0("max weight: ", round(max(Weight), 2)))
  } else {
    compare_distrib <- data.table("Coupling angle rescale" = scales::rescale(edges_angle$weight, to = c(0,1)),
                                  "Coupling strength rescale" = scales::rescale(edges_strength$weight, to = c(0,1)),
                                  "Coupling similarity rescale" = scales::rescale(edges_similarity$weight, to = c(0,1)))
    compare_distrib <- pivot_longer(compare_distrib, cols = 1:3, names_to = "Method", values_to = "Weight") %>% 
      group_by(Method) %>% 
      mutate(label = paste0("max weight: ", round(max(Weight), 2)))
  }
  
  plot <- ggplot(compare_distrib, aes(Weight, after_stat(count), fill = Method)) +
    geom_density(show.legend = FALSE) +
    geom_text(aes(label = label, x = Inf, y = Inf), hjust = 1.5, vjust = 1.5) +
    facet_wrap(~ Method, scales = "free", ncol = 1) +
    theme_light()
}

compare_partition <- function(graph, resolutions = seq(0.5,1.5, 0.1)){
  #' Compare Partitions for Different Leiden Algorithm Resolutions
  #'
  #' This function runs the Leiden algorithm for different resolutions value and
  #' compares the distribution of communities regarding their size, the number 
  #' of communities and the modularity
  #'
  #' @param graph
  #' A tidygraph network.
  #'
  #' @param resolutions
  #' Choose the different resolutions you want to test with the Leiden algorithm.
  #' By default, the function starts from 0.5 and goes to 1.5 by decile.   
  
  stat_partition <- data.table("Com_ID" = c(),
                               "Size_Com" = c(),
                               "Modularity" = c(),
                               "Resolution" = c(),
                               "n_community" = c())
  for(i in resolutions){
    graph <- networkflow::leiden_workflow(graph, res_1 = i, niter = 1000)
    
    partition <- graph %>% 
      activate(nodes) %>% 
      as.data.table() %>% 
      select(Com_ID, Size_com) %>%
      unique() %>% 
      mutate(Modularity = modularity(graph, membership = V(graph)$Com_ID),
             Resolution = i,
             n_community = length(unique(V(graph)$Com_ID)))
    
    stat_partition <- rbind(stat_partition, partition)
  }
  
  coeff_y_axis <- max(stat_partition$n_community) * 1.1
  
  stat_partition %>% 
    ggplot(aes(x = as.character(Resolution), y = Size_com, fill = Com_ID, group = Com_ID)) +
    geom_col(position = "fill") +
    scale_fill_viridis_d(option = "magma", direction = -1) +
    geom_point(aes(y = n_community / coeff_y_axis), size = 2, show.legend = FALSE) +
    geom_line(aes(y = Modularity), size = 1.5) +
    scale_y_continuous(name = "Community size and modularity", 
                       sec.axis = sec_axis(~.*coeff_y_axis, name = "Number of communities")) +
    scale_x_discrete(name = "Resolutions") +
    guides(fill = guide_legend(title = "Communities")) +
    geom_vline(aes(xintercept = as.character(unique(stat_partition[Modularity == max(Modularity)]$Resolution))), linetype = "dashed") +
    theme_minimal()
}

name_community_with_words <- function(graph, title_column = "Title", com_column = "Com_ID", n_words = 4, remove_words = NULL){
  #' Automatically Attributing Names to Communities with Top TF-IDF Words
  #'
  #' @description A function to give to a community the name of the most identifying words
  #' in titles and/or abstracts of the documents within the community.
  #'
  #' @param graph A tidygraph object.
  #'
  #' @param title_column Enter the name of the column with the words that will
  #' be used for naming the communities.
  #'
  #' @param com_column
  #' The name of your community identifiers column.
  #' 
  #' @param n_words
  #' The number of words you want to have in the name of the community.
  #'
  #' @details The attribute of nodes and edges with the names of the communities is called
  #' `com_name_by_words`.
  #'
  #' @return The same graph object but with a column `com_name_by_words`.

# Finding the words with the highest tf-idf value per community
  # extracting the nodes
    nodes <- graph %>%
      activate(nodes) %>%
      as.data.table()

    # changing the names of the column for titles and communities
    colnames(nodes)[colnames(nodes) == com_column] <- "Com_ID"
    colnames(nodes)[colnames(nodes) == title_column] <- "Titre"
    
    tf_idf <- nodes[Titre != "NULL"] %>% 
      select(Com_ID, Titre) %>% 
      tidytext::unnest_tokens(word, Titre) %>% 
      anti_join(tidytext::stop_words) %>% 
      mutate(word := textstem::lemmatize_words(word)) %>% 
      count(Com_ID, word) %>%
      tidytext::bind_tf_idf(word, Com_ID, n)
    
    if(length(remove_words) > 0){
    tf_idf <- tf_idf %>% 
      filter(! word %in% remove_words)
    } 
    
    top_words <- tf_idf %>% 
      group_by(Com_ID) %>% 
      slice_max(order_by = tf_idf, n = n_words, with_ties = FALSE) %>% 
      mutate(com_name_by_words = paste0(word, collapse = "-")) %>% 
      select(Com_ID, com_name_by_words) %>% 
      unique()
    
    # adding the name as an attribute to the nodes.
    graph_coupling <- graph_coupling %>%
      activate(nodes) %>%
      inner_join(top_words, by = "Com_ID")
}

label_com <- function(graph, biggest_community = FALSE, community_threshold = 0.01, community_name_column = "Community_name", community_size_column = "Size_com") {
  #' Displaying the highest cited nodes
  #'
  #' A simple function to calculate the mean of coordinates x and y for each community. These coordinates
  #' are to be used to plot the name of the community on the graph.
  #'
  #' @param graph
  #' A tidygraph object.
  #'
  #' @param biggest_community
  #' If true, you have the possibility to remove the smallest community, depending of the `community_threshold`
  #' you have set.
  #'
  #' @param community_threshold
  #' If `biggest_community` is true, the function selects the nodes that belong to communities which represent
  #' at least x% of the total number of nodes. By default, the parameter is set to 1%.
  #'
  #' @param community_name_column
  #' Name of the column with the name of the community to be used as label
  #'
  #' @param community_size_column
  #' Name of the column with the share of nodes in each community.
  
  
  # Top nodes per community for the variable chosen
  label_com <- graph %>%
    activate(nodes) %>%
    as_tibble()
  
  # Changing the name of the variable chosen
  colnames(label_com)[colnames(label_com) == community_name_column] <- "Community_name"
  colnames(label_com)[colnames(label_com) == community_size_column] <- "Size_com"
  
  # Keeping only the biggest communites if the parameter is TRUE
  if (biggest_community == TRUE) {
    label_com <- label_com %>%
      filter(Size_com > community_threshold)
  }
  
  # Keeping the n top nodes per community
  label_com <- label_com %>%
    group_by(Community_name) %>%
    mutate(x = mean(x), y = mean(y)) %>%
    select(Community_name, x, y, color, Size_com) %>%
    as_tibble() %>%
    unique()
  
  return(label_com)
}

#' # Building Graphs - Secondary and to be improved

#' This function takes as input a tidygraph object, it switches the names of individual nodes
#' by the name of their community, calculate the number of nodes in the community, and
#' transforms the whole community as a unique node.
#' 
graph_community <- function(graph, Community_name = "Community_name", nb_components = 1, preparing_graph = TRUE) {
  #' Function for building of graph with community as nodes
  #'
  #' This function takes as input a tidygraph object, it switches the names of individual nodes
  #' by the name of their community, calculate the number of nodes in the community, and
  #' transforms the whole community as a unique node.
  #'
  #' @param graph
  #' A tidygraph object.
  #'
  #' @param Community_name
  #' Select the column that you want to be used to name the community. By default, the function
  #' considers that you have a column called `Community_name`. Be careful to choose a variable that
  #' identifies all members of a community share. It could be the number of the community, or the
  #' name of the node with the highest degree, number of citations, etc...
  #'
  #' @param nb_components
  #' The number of components you want to keep in the network. Usually, as nodes are community, it is
  #' likely that no node is isolated from the rest of the network.
  #'
  #' @param preparing_graph
  #' If `TRUE`, the function prepares the graph to be plotted with ggraph. It attributes color to the name
  #' of the nodes, depending of the color attributed in the original graph (the tidygraph use in the parameter
  #' `graph`). The `graph` which serves as input thus need to already have a `color` column. The function
  #' then run the Force Atlas algorithm to calculate the position of nodes.
  #'
  #' @details
  #' For running the Force Atlas, we have standard parameters that are not modifiable via the parameters of the
  #' function, as in general one doesn't have more than 40/50 communities and that you don't need
  #' to adjust the parameters. If one has many communities and one wants to find the communities of these
  #' communities as node, it is possible to do out from the output of the function.
  #'
  #' @details
  #' The links between nodes(ie communities) are normalized such as not being overestimated for bigger
  # communities with more edges.
  
  
  Nodes_Com <- graph %>%
    activate(nodes) %>%
    as_tibble()
  
  colnames(Nodes_Com)[colnames(Nodes_Com) == Community_name] <- "Community_name"
  
  Nodes_Com <- Nodes_Com %>%
    mutate(name = Community_name) %>%
    select(name) %>%
    group_by(name) %>%
    mutate(nb_nodes = n()) %>%
    unique()
  
  # associating the two nodes of each edge with their respective community name
  # summing the weights of the similar edges
  Edges_Com <- graph %>%
    activate(edges) %>%
    mutate(Com_name_from = .N()$Community_name[from], Com_name_to = .N()$Community_name[to]) %>%
    as_tibble() %>%
    select(Com_name_from, Com_name_to, weight) %>%
    group_by(Com_name_from, Com_name_to) %>%
    mutate(weight = sum(weight)) %>%
    unique() %>%
    as.data.table()
  
  # calculating the sum of weigths for each node (i.e. each community), that is the sum
  # of the weights of all the node edges
  Weight_com1 <- Edges_Com[Com_name_from != Com_name_to, c("Com_name_from", "weight")]
  Weight_com2 <- Edges_Com[Com_name_from != Com_name_to, c("Com_name_to", "weight")]
  Weight_com3 <- Edges_Com[Com_name_from == Com_name_to, c("Com_name_from", "weight")]
  Weight_com <- rbind(Weight_com1, Weight_com2, Weight_com3, use.names = FALSE)
  Weight_com <- Weight_com[, Com_weight := sum(weight), by = "Com_name_from"]
  colnames(Weight_com)[1] <- "Com_name"
  Weight_com <- unique(Weight_com[, c("Com_name", "Com_weight")])
  
  # merging the two nodes of each edge with their respective total weight
  Edges_Com <- merge(Edges_Com, Weight_com, by.x = "Com_name_from", by.y = "Com_name")
  Edges_Com <- merge(Edges_Com, Weight_com, by.x = "Com_name_to", by.y = "Com_name")
  Edges_Com <- Edges_Com[, c("Com_name_from", "Com_name_to", "weight", "Com_weight.x", "Com_weight.y")]
  colnames(Edges_Com) <- c("from", "to", "weight", "Com_weight_from", "Com_weight_to")
  
  # Cosine normalization of the edges weights
  Edges_Com <- Edges_Com[, weight := weight / sqrt(Com_weight_from * Com_weight_to), by = c("from", "to")]
  
  # using our tbl_main_components function to build the tidygraph with a main component
  graph_community <- tbl_main_components(Edges_Com, Nodes_Com, node_key = "name", nb_components = nb_components)
  
  if (preparing_graph == TRUE) {
    # merging with the colors attributed to community before
    color_community <- graph %>%
      activate(nodes) %>%
      as_tibble() %>%
      select(Community_name, color) %>%
      unique()
    
    graph_community <- graph_community %>%
      activate(nodes) %>%
      left_join(color_community, by = c("Id" = "Community_name"))
    
    
    # Integration a size variable for implementing non-overlapping function of Force Atlas
    graph_community <- graph_community %>%
      activate(nodes) %>%
      mutate(size = nb_nodes)
    
    # Running Force Atlas layout
    graph_community <- force_atlas(graph_community, seed = NULL, ew.influence = 1, kgrav = 1, iter_1 = 6000, iter_2 = 2000, barnes.hut = FALSE, size_min = 50, size_max = 200)
  }
}


# A function to concentrate nodes with a similar attribute in one singular node and build the corresponding graph
graph_from_attribute <- function(nodes, edges, palette, Attribute_name, nb_components = 1, preparing_graph = TRUE, size_min = 50, size_max = 200) {
  #' Function for building of graph with an attribute as nodes
  #'
  #' This function takes as input a tidygraph object, it switches the names of individual nodes
  #' by one of their attribute, calculate the number of nodes having this attribute, and
  #' transforms all the nodes with this attribute as a unique node.
  #'
  #' @param nodes
  #' The list of the nodes of your network.
  #'
  #' @param edges
  #' The list of the edges of your network.
  #'
  #' @param palette
  #' A color palette that will be used to attribute color to the communities calculated
  #' by the Leiden algorithm
  #'
  #' @param Attribute_name
  #' Select the column that you want to be used for aggregating the nodes.
  #'
  #' @param nb_components
  #' The number of components you want to keep in the network. Usually, as nodes are community, it is
  #' likely that no node is isolated from the rest of the network.
  #'
  #' @param preparing_graph
  #' If `TRUE`, the function prepares the graph to be plotted with ggraph. It will run Leiden, attribute
  #' color to the community, calculate some centrality measures and then run the Force Atlast
  #' algorithm to calculate the position of nodes.
  #'
  #' @param size_min
  #' A parameter used in the Force Atlas algorithm to avoid the overlappping of nodes. See the
  #' `force_atlas` function for documentation.
  #'
  #' @param size_max
  #' See `size_min`.
  #'
  
  colnames(nodes)[colnames(nodes) == Attribute_name] <- "Attribute_name"
  graph <- tbl_main_components(edges = edges, nodes = nodes, node_key = "ID_Art", threshold_alert = 0.05, directed = FALSE)
  
  Nodes_Com <- graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(name = Attribute_name) %>%
    select(name) %>%
    group_by(name) %>%
    mutate(Size_att = n()) %>%
    arrange(desc(Size_att)) %>%
    unique()
  
  # associating the two nodes of each edge with their respective community name
  # summing the weights of the similar edges
  Edges_Com <- graph %>%
    activate(edges) %>%
    mutate(Att_name_from = .N()$Attribute_name[from], Att_name_to = .N()$Attribute_name[to]) %>%
    as_tibble() %>%
    select(Att_name_from, Att_name_to, weight) %>%
    group_by(Att_name_from, Att_name_to) %>%
    mutate(weight = sum(weight)) %>%
    unique() %>%
    as.data.table()
  
  # calculating the sum of weigths for each node (i.e. each community), that is the sum
  # of the weights of all the node edges
  Weight_com1 <- Edges_Com[Att_name_from != Att_name_to, c("Att_name_from", "weight")]
  Weight_com2 <- Edges_Com[Att_name_from != Att_name_to, c("Att_name_to", "weight")]
  Weight_com3 <- Edges_Com[Att_name_from == Att_name_to, c("Att_name_from", "weight")]
  Weight_com <- rbind(Weight_com1, Weight_com2, Weight_com3, use.names = FALSE)
  Weight_com <- Weight_com[, Com_weight := sum(weight), by = "Att_name_from"]
  colnames(Weight_com)[1] <- "Att_name"
  Weight_com <- unique(Weight_com[, c("Att_name", "Com_weight")])
  
  # merging the two nodes of each edge with their respective total weight
  Edges_Com <- merge(Edges_Com, Weight_com, by.x = "Att_name_from", by.y = "Att_name")
  Edges_Com <- merge(Edges_Com, Weight_com, by.x = "Att_name_to", by.y = "Att_name")
  Edges_Com <- Edges_Com[, c("Att_name_from", "Att_name_to", "weight", "Com_weight.x", "Com_weight.y")]
  colnames(Edges_Com) <- c("from", "to", "weight", "Att_name_from", "Att_name_to")
  
  # Cosine normalization of the edges weights
  Edges_Com <- Edges_Com[, weight := weight / sqrt(Att_name_from * Att_name_to), by = c("from", "to")]
  
  # using our tbl_main_components function to build the tidygraph with a main component
  graph_community <- tbl_main_components(Edges_Com, Nodes_Com, node_key = "name", nb_components = nb_components)
  
  if (preparing_graph == TRUE) {
    # Identifying communities with Leiden algorithm
    graph_community <- leiden_improved(graph_community, res_1 = 1, res_2 = NULL, res_3 = NULL, n_iterations = 500)
    
    # Giving colors to communities
    graph_community <- community_colors(graph_community, palette)
    
    # Calculating different centrality measures
    graph_community <- centrality(graph_community)
    
    # Integration a size variable for implementing non-overlapping function of Force Atlas
    graph_community <- graph_community %>%
      activate(nodes) %>%
      mutate(size = Size_att)
    
    # Running Force Atlas layout
    graph_community <- force_atlas(graph_community, seed = NULL, ew.influence = 1, kgrav = 1, iter_1 = 6000, iter_2 = 2000, barnes.hut = TRUE, size_min = size_min, size_max = size_max)
  }
}

clustering_communities <- function(graph, label_size = 6, number_size = 6, threshold_com = 0.01) {
  #' Function for building a heatmap of the communities
  #'
  #' This function takes as input a tidygraph object with communities as nodes and produce a heatmap of the links
  #' between communities, and a dendrogram of these communities.
  #'
  #' @param graph
  #' A tidygraph object with nodes being communities and a column `Size_com` which represent the percentage of
  #' total nodes in the community
  #'
  #' @param label_size
  #' The size of the labels displayed in the heatmap plot
  #'
  #' #' @param threshold_com
  #' The minimun percentage of nodes in the community for the community to be displayed on the plot.
  #'
  #' @section Future improvements:
  #' Find a way to plot only the biggest communities, but by removing the communities before the plotting, not to
  #' biased the values.
  #'
  #' @section Future improvements:
  #' Using a different method for plotting, by mixing ggplot and ggraph for the dendrogram, to have more options.
  
  # Extracting edges with only the biggest communities, and integrating the name of communities for source and target of edges
  edges <- graph %>%
    activate(edges) %>%
    mutate(com_name_to = .N()$Id[to], com_name_from = .N()$Id[from]) %>%
    as.data.table()
  
  # making the matrix from the edges
  matrix <- as.matrix(get.adjacency(graph.data.frame(edges[, c("com_name_to", "com_name_from", "weight")], directed = FALSE), type = "both", attr = "weight"))
  
  # clustering and creation of a dendrogram from the Matrix.
  dendro <- as.dendrogram(hclust(dist(matrix)))
  plot_dendro <- ggdendrogram(dendro, rotate = TRUE) # saving the plot of the dendrogram
  order_dendro <- order.dendrogram(dendro) # extracting the order of nodes
  
  # keeping only biggest communities
  nodes <- graph %>%
    activate(nodes) %>%
    as.data.table()
  
  nodes <- nodes[, total := sum(nb_nodes)][, share_com := nb_nodes / total][share_com > threshold_com, "Id"]
  
  # cleaning the matrix for plotting the heatmap
  matrix <- scale(matrix, center = FALSE, scale = colSums(matrix))
  matrix <- melt(matrix) %>% as.data.table()
  matrix$value <- matrix$value * 100
  matrix$value <- round(matrix$value, digits = 1)
  matrix <- matrix[matrix$Var1 %in% nodes$Id & matrix$Var2 %in% nodes$Id]
  matrix$Var1 <- str_wrap(matrix$Var1, width = 10)
  matrix$Var2 <- str_wrap(matrix$Var2, width = 200)
  
  # ordering the nodes depending of the order of the dendrogram
  matrix$Var1 <- factor(
    x = matrix$Var1,
    levels = unique(matrix$Var1)[order_dendro],
    ordered = TRUE
  )
  matrix$Var2 <- factor(
    x = matrix$Var2,
    levels = unique(matrix$Var2)[order_dendro],
    ordered = TRUE
  )
  
  # saving the heat map
  plot_heatmap <- ggplot(matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(show.legend = FALSE) +
    theme(text = element_text(size = label_size)) +
    geom_text(aes(x = Var1, y = Var2, label = value), color = "black", size = number_size) +
    scale_fill_viridis(discrete = FALSE) +
    ylab("In the cluster...") +
    xlab("...X% of links goes to")
  
  list_return <- list("heatmap" = plot_heatmap, "dendrogram" = plot_dendro)
  return(list_return)
}


#' # Functions for word analysis (titles) of networks 
#' 
#' ## `tf_idf()`
#' 
#' This function takes as input a tidygraph object or a data frame with nodes, both with a community attribute, and analyzes
#' the words use in the title of the articles to calculate the words with the highest TF-IDF
#' value for each community.
#' 

tf_idf <- function(graph = NULL, nodes = NULL, title_column = "Titre", com_column = "Com_ID", color_column = "color",
                   com_name_column = "Community_name", com_size_column = "Size_com", threshold_com = 0.01, number_of_words = 12,
                   palette = NULL, size_title_wrap = 8, lemmatize_bigrams = TRUE, plot = TRUE) {
  #' Creating a TF-IDF analysis of the titles of WoS corpus
  #'
  #' This function takes as input a tidygraph object or a data frame with nodes, both with a community attribute, and analyzes
  #' the words use in the title of the articles to calculate the words with the highest TF-IDF
  #' value for each community.
  #'
  #' @param graph
  #' A tidygraph object. By default `NULL` in case you prefer to enter a data frame with nodes.
  #'
  #' @param nodes
  #' A data frame with the nodes of the network, and community and title attributes.
  #'
  #' @param title_column
  #' The name of the column with the titles of the articles. The function renames the column
  #' "Titre", as in the OST WoS database ("Titre" is the default value).
  #'
  #' @param com_column
  #' The name of the column with the id of the communities. The function renames the column
  #' "Com_ID" (default value).
  #'
  #' @param color_column
  #' The name of the column with the color attribute of the communities. The function renames the column
  #' "color" (default value).
  #'
  #' @param com_name_column
  #' The name of the column with the name of the communities.
  #'
  #' @param com_size_column
  #' The name of the column with the share of total nodes in each community.
  #'
  #' @param threshold_com
  #' The minimun percentage of nodes in the community for the community to be displayed on the plot.
  #'
  #' @param number_of_words
  #' How many words you want to display on the final graph.
  #'
  #' @param palette
  #' If you don't already have a color attribute for your communities in your tidygraph object,
  #' the function will generate one from a palette that you can add in the paramaters (NULL by default).
  #'
  #' @param size_title_wrap
  #' The size of the community title in the plot.
  #' 
  #' @param lemmatize_bigrams
  #' Chose whether you want to lemmatize each word in a bigram (which could lead to a bigram that has no 
  #' clear meaning) or not.
  
  # extracting the nodes
  if (!is.null(graph)) {
    tf_idf_save <- graph %>%
      activate(nodes) %>%
      as.data.table()
  }
  else {
    tf_idf_save <- nodes %>% as.data.table()
  }
  
  # changing the names of the column for titles and communities
  colnames(tf_idf_save)[colnames(tf_idf_save) == com_column] <- "Com_ID"
  colnames(tf_idf_save)[colnames(tf_idf_save) == title_column] <- "Titre"
  colnames(tf_idf_save)[colnames(tf_idf_save) == com_name_column] <- "Community_name"
  colnames(tf_idf_save)[colnames(tf_idf_save) == com_size_column] <- "Size_com"
  
  
  # adding a color column attribute in case it doesn't exist
  if (!any(names(tf_idf_save) == color_column)) {
    color <- data.table(
      Com_ID = 1:500,
      color = palette
    )
    color <- color %>%
      mutate(Com_ID = sprintf("%02d", Com_ID)) %>%
      mutate(Com_ID = as.character(Com_ID))
    
    tf_idf <- merge(tf_idf_save, color, by = "Com_ID", all.x = TRUE)
  } else {
    colnames(tf_idf_save)[colnames(tf_idf_save) == color_column] <- "color"
  }
  
  tf_idf <- tf_idf_save # we will need tf_idf_save later
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  #### Unigram ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # Cleaning the titles
  
  tf_idf <- tf_idf[Titre != "NULL"]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := gsub("-", " ", Titre)] # to avoid having non-existing words like UnitedStates
  tf_idf[, Titre := stripWhitespace(Titre)]
  tf_idf[, Titre := removePunctuation(Titre)]
  tf_idf[, Titre := removeNumbers(Titre)]
  tf_idf[, Titre := tolower(Titre)]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := as.character(Titre)]
  tf_idf$Titre <- quanteda::tokens(tf_idf$Titre, remove_punct = TRUE)
  tf_idf$Titre <- quanteda::tokens_ngrams(tf_idf$Titre, n = 1)
  tible_tf_idf <- tf_idf[, paste(Titre, collapse = " "), by = "Com_ID"]
  tible_tf_idf[, V1 := stripWhitespace(V1)]
  # Dictionnary to find the root of stem word before stemming
  tf_idf_table <- tible_tf_idf
  tf_idf_table <- tf_idf_table %>% unnest_tokens(word, V1) %>% as.data.table()
  tf_idf_table <- tf_idf_table[, word := textstem::lemmatize_words(word)] 
  tf_idf_table <- tf_idf_table[, count := .N, by = c("Com_ID","word")] %>% unique()
  # applying tf-idf
  tf_idf_table <- tidytext::bind_tf_idf(tf_idf_table, word, Com_ID, count)
  tf_idf_table_uni <- tf_idf_table
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  #### Bigram ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # cleaning
  
  tf_idf <- tf_idf_save[Titre != "NULL"]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := gsub("-", " ", Titre)] # to avoid having non-existing words like UnitedStates
  tf_idf[, Titre := stripWhitespace(Titre)]
  tf_idf[, Titre := removePunctuation(Titre)]
  tf_idf[, Titre := removeNumbers(Titre)]
  tf_idf[, Titre := tolower(Titre)]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := as.character(Titre)]
  # ngraming
  # tf_idf[, Titre := stemDocument(Titre)]
  tf_idf$Titre <- tokens(tf_idf$Titre, remove_punct = TRUE)
  tf_idf$Titre <- tokens_ngrams(tf_idf$Titre, n = 2)
  tible_tf_idf <- tf_idf[, paste(Titre, collapse = " "), by = "Com_ID"]
  tible_tf_idf[, V1 := stripWhitespace(V1)]
  
  tf_idf_table <- tible_tf_idf
  tf_idf_table <- tf_idf_table %>% unnest_tokens(word, V1) %>% as.data.table()
  tf_idf_table$word <- gsub("_", " ", tf_idf_table$word)
  
  if(lemmatize_bigrams == TRUE){
    tf_idf_table[, word1 := str_extract(tf_idf_table$word, "\\S+")]
    tf_idf_table[, word2 := str_extract(tf_idf_table$word, "\\S+$")]
    tf_idf_table <- tf_idf_table[, word1 := textstem::lemmatize_words(word1)]
    tf_idf_table <- tf_idf_table[, word2 := textstem::lemmatize_words(word2)] 
    tf_idf_table <- tf_idf_table %>%
      tidyr::unite(word, word1, word2, sep = " ") %>% 
      as.data.table()
  }
  
  tf_idf_table <- tf_idf_table[, count := .N, by = c("Com_ID","word")] %>% unique()
  
  # applying tf-idf
  tf_idf_table <- tidytext::bind_tf_idf(tf_idf_table, word, Com_ID, count)
  tf_idf_table_bi <- tf_idf_table
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  #### Plot ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  tf_idf_table <- rbind(tf_idf_table_uni, tf_idf_table_bi, fill = TRUE)
  tf_idf_table <- tf_idf_table[order(-tf_idf)][, head(.SD, number_of_words), Com_ID]
  
  # Get info about size of communities
  Size_com <- unique(tf_idf_save[, .(Com_ID, Community_name, Size_com, color)])
  tf_idf_table <- merge(tf_idf_table, Size_com, by = "Com_ID", all.x = TRUE) # merge
  
  # Wrap Name and Reorder according to share_leiden
  tf_idf_table$Com_wrap <- str_wrap(tf_idf_table$Community_name, width = 10)
  tf_idf_table <- tf_idf_table[order(-Size_com)] # order by variable
  tf_idf_table <- tf_idf_table[order(-Size_com)] # order by variable
  tf_idf_table$Com_wrap <- factor(tf_idf_table$Com_wrap) # make a factor
  tf_idf_table$Com_wrap <- fct_inorder(tf_idf_table$Com_wrap) # by order of appearance
  
  # fixing the number of columns, depending of the number of communities
  n_columns <- 3
  if (length(unique(tf_idf_table[Size_com >= threshold_com]$Com_wrap)) > 9) {
    n_columns <- 4
  }
  
  if (length(unique(tf_idf_table[Size_com >= threshold_com]$Com_wrap)) > 12) {
    n_columns <- 5
  }
  
  if (length(unique(tf_idf_table[Size_com >= threshold_com]$Com_wrap)) > 15) {
    n_columns <- 6
  }
  
  # plotting the graph
  if (plot == TRUE){
  tf_idf_plot <- ggplot(tf_idf_table[Size_com >= threshold_com], aes(reorder_within(word, tf_idf, color), tf_idf, fill = color)) +
    geom_bar(stat = "identity", alpha = .8, show.legend = FALSE) +
    labs(
      title = "Highest tf-idf",
      x = "Words", y = "tf-idf"
    ) +
    facet_wrap(~Com_wrap, ncol = n_columns, scales = "free") +
    scale_x_reordered() +
    scale_fill_identity() +
    theme(strip.text = element_text(size = size_title_wrap)) +
    coord_flip()
  
  list_return <- list("plot" = tf_idf_plot, "list_words" = tf_idf_table)
  return(list_return)
  } else {
    return(tf_idf_table)
  }
}


#' ## `tf_idf_old()`
#' 
#' > This is former function that used stemming. We now prefer using lemmatization.
#' 
#' 
#' This function takes as input a tidygraph object or a data frame with nodes, both with a community attribute, and analyzes
#' the words use in the title of the articles to calculate the words with the highest TF-IDF
#' value for each community.
#' 

tf_idf_stemmed <- function(graph = NULL, nodes = NULL, title_column = "Titre", com_column = "Com_ID", color_column = "color",
                           com_name_column = "Community_name", com_size_column = "Size_com", threshold_com = 0.01, number_of_words = 12,
                           palette = NULL, size_title_wrap = 8, unstemming = TRUE) {
  #' Creating a TF-IDF analysis of the titles of WoS corpus
  #'
  #' This function takes as input a tidygraph object or a data frame with nodes, both with a community attribute, and analyzes
  #' the words use in the title of the articles to calculate the words with the highest TF-IDF
  #' value for each community.
  #'
  #' @param graph
  #' A tidygraph object. By default `NULL` in case you prefer to enter a data frame with nodes.
  #'
  #' @param nodes
  #' A data frame with the nodes of the network, and community and title attributes.
  #'
  #' @param title_column
  #' The name of the column with the titles of the articles. The function renames the column
  #' "Titre", as in the OST WoS database ("Titre" is the default value).
  #'
  #' @param com_column
  #' The name of the column with the id of the communities. The function renames the column
  #' "Com_ID" (default value).
  #'
  #' @param color_column
  #' The name of the column with the color attribute of the communities. The function renames the column
  #' "color" (default value).
  #'
  #' @param com_name_column
  #' The name of the column with the name of the communities.
  #'
  #' @param com_size_column
  #' The name of the column with the share of total nodes in each community.
  #'
  #' @param threshold_com
  #' The minimun percentage of nodes in the community for the community to be displayed on the plot.
  #'
  #' @param number_of_words
  #' How many words you want to display on the final graph.
  #'
  #' @param palette
  #' If you don't already have a color attribute for your communities in your tidygraph object,
  #' the function will generate one from a palette that you can add in the paramaters (NULL by default).
  #'
  #' @param size_title_wrap
  #' The size of the community title in the plot.
  #' 
  #' @param unstemming
  #' Chose whether you want to unstem or not the words
  
  # extracting the nodes
  if (!is.null(graph)) {
    tf_idf_save <- graph %>%
      activate(nodes) %>%
      as.data.table()
  }
  else {
    tf_idf_save <- nodes %>% as.data.table()
  }
  
  # changing the names of the column for titles and communities
  colnames(tf_idf_save)[colnames(tf_idf_save) == com_column] <- "Com_ID"
  colnames(tf_idf_save)[colnames(tf_idf_save) == title_column] <- "Titre"
  colnames(tf_idf_save)[colnames(tf_idf_save) == com_name_column] <- "Community_name"
  colnames(tf_idf_save)[colnames(tf_idf_save) == com_size_column] <- "Size_com"
  
  
  
  # adding a color column attribute in case it doesn't exist
  if (any(names(tf_idf_save) == color_column)) {
    color <- data.table(
      Com_ID = 1:500,
      color = mypalette
    )
    color <- color %>%
      mutate(Com_ID = sprintf("%02d", Com_ID)) %>%
      mutate(Com_ID = as.character(Com_ID))
    
    tf_idf <- merge(tf_idf_save, color, by = "Com_ID", all.x = TRUE)
  } else {
    colnames(tf_idf_save)[colnames(tf_idf_save) == color_column] <- "color"
  }
  
  tf_idf <- tf_idf_save # we will need tf_idf_save later
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  #### Unigram ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # Cleaning the titles
  
  tf_idf <- tf_idf[Titre != "NULL"]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := stripWhitespace(Titre)]
  tf_idf[, Titre := removePunctuation(Titre)]
  tf_idf[, Titre := removeNumbers(Titre)]
  tf_idf[, Titre := tolower(Titre)]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := as.character(Titre)]
  tf_idf$Titre <- tokens(tf_idf$Titre, remove_punct = TRUE)
  tf_idf$Titre <- tokens_ngrams(tf_idf$Titre, n = 1)
  tible_tf_idf <- tf_idf[, paste(Titre, collapse = " "), by = "Com_ID"]
  tible_tf_idf[, V1 := stripWhitespace(V1)]
  # Dictionnary to find the root of stem word before stemming
  dictionary <- tible_tf_idf
  dictionary <- dictionary %>% unnest_tokens(word, V1) %>% as.data.table()
  tible_tf_idf[, V1 := stemDocument(V1)]
  # tf-idf using quanteda
  tible_tf_idf <- corpus(tible_tf_idf, text_field = "V1")
  tible_tf_idf <- dfm(tible_tf_idf)
  tible_tf_idf <- dfm_tfidf(tible_tf_idf)
  # Keep column of names
  documents_names <- cbind(docvars(tible_tf_idf), quanteda::convert(tible_tf_idf, to = "data.frame")) %>% as.data.table()
  tible_tf_idf <- tidy(tible_tf_idf) %>% as.data.table()
  tf_idf_table <- merge(tible_tf_idf, documents_names[, .(doc_id, Com_ID)], by.x = "document", by.y = "doc_id")
  setkey(tf_idf_table,term)
  setkey(dictionary, word)
  
  if(unstemming==TRUE){
    # we separate the terms to optimise the matching with unstemmed words
    terms <- unique(tf_idf_table$term) %>% as.data.table()
    setnames(terms,".","term")
    setkey(terms, term)
    terms <- terms[, unstemmed_word:= stemCompletion(term, dictionary$word, type = "prevalent")]
    tf_idf_table <- merge(tf_idf_table, terms, by = "term")
    tf_idf_table[unstemmed_word=="",unstemmed_word:=term] # unstem with most common word
  }
  if(unstemming==FALSE){
    tf_idf_table[,unstemmed_word:=term]
  }
  
  tf_idf_table[, term := unstemmed_word]
  tf_idf_table_uni <- tf_idf_table
  
  
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  #### Bigram ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # cleaning
  
  tf_idf <- tf_idf_save[Titre != "NULL"]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := stripWhitespace(Titre)]
  tf_idf[, Titre := removePunctuation(Titre)]
  tf_idf[, Titre := removeNumbers(Titre)]
  tf_idf[, Titre := tolower(Titre)]
  tf_idf[, Titre := removeWords(Titre, stopwords("english"))]
  tf_idf[, Titre := as.character(Titre)]
  # ngraming
  tf_idf[, Titre := stemDocument(Titre)]
  tf_idf$Titre <- tokens(tf_idf$Titre, remove_punct = TRUE)
  tf_idf$Titre <- tokens_ngrams(tf_idf$Titre, n = 2)
  tible_tf_idf <- tf_idf[, paste(Titre, collapse = " "), by = "Com_ID"]
  tible_tf_idf[, V1 := stripWhitespace(V1)]
  # tf-idf using quanteda
  tible_tf_idf <- corpus(tible_tf_idf, text_field = "V1")
  tible_tf_idf <- dfm(tible_tf_idf)
  tible_tf_idf <- dfm_tfidf(tible_tf_idf)
  # Keep column of names
  documents_names <- cbind(docvars(tible_tf_idf), quanteda::convert(tible_tf_idf, to = "data.frame")) %>% as.data.table()
  tible_tf_idf <- tidy(tible_tf_idf) %>% as.data.table()
  tf_idf_table <- merge(tible_tf_idf, documents_names[, .(doc_id, Com_ID)], by.x = "document", by.y = "doc_id")
  # Unstemming bigram: first term, then second term, them bringing them together
  tf_idf_table$term <- gsub("_", " ", tf_idf_table$term)
  tf_idf_table[, term1 := str_extract(tf_idf_table$term, "\\S+")]
  tf_idf_table[, term2 := str_extract(tf_idf_table$term, "\\S+$")]
  setkey(tf_idf_table,term1, term2)
  setkey(dictionary, word)
  
  if(unstemming==TRUE){
    terms <- unique(tf_idf_table$term1) %>% as.data.table()
    terms2 <- unique(tf_idf_table$term2) %>% as.data.table()
    terms <- rbind(terms,terms2) %>% unique()
    setnames(terms,".","term")
    setkey(terms, term)
    terms <- terms[, unstemmed_word:= stemCompletion(term, dictionary$word, type = "prevalent")]
    tf_idf_table <- merge(tf_idf_table, terms, by.x = "term1", by.y = "term")
    tf_idf_table <- merge(tf_idf_table, terms, by.x = "term2", by.y = "term")
    tf_idf_table[unstemmed_word.x=="",unstemmed_word:=term1]
    tf_idf_table[unstemmed_word.y=="",unstemmed_word:=term2]
    setnames(tf_idf_table,c("unstemmed_word.x","unstemmed_word.y"),c("unstemmed_word1","unstemmed_word2"))
  }
  if(unstemming==FALSE){
    tf_idf_table[,unstemmed_word1:=term1]
    tf_idf_table[,unstemmed_word2:=term2]
  }
  
  tf_idf_table[, term := paste(unstemmed_word1, unstemmed_word2)]
  tf_idf_table_bi <- tf_idf_table
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  #### Plot ####
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  tf_idf_table <- rbind(tf_idf_table_uni, tf_idf_table_bi, fill = TRUE)
  tf_idf_table <- tf_idf_table[order(-count)][, head(.SD, number_of_words), Com_ID]
  
  # Get info about size of communities
  Size_com <- unique(tf_idf_save[, .(Com_ID, Community_name, Size_com, color)])
  tf_idf_table <- merge(tf_idf_table, Size_com, by = "Com_ID", all.x = TRUE) # merge
  
  # Wrap Name and Reorder according to share_leiden
  tf_idf_table$Com_wrap <- str_wrap(tf_idf_table$Community_name, width = 10)
  tf_idf_table <- tf_idf_table[order(-Size_com)] # order by variable
  tf_idf_table <- tf_idf_table[order(-Size_com)] # order by variable
  tf_idf_table$Com_wrap <- factor(tf_idf_table$Com_wrap) # make a factor
  tf_idf_table$Com_wrap <- fct_inorder(tf_idf_table$Com_wrap) # by order of appearance
  
  # fixing the number of columns, depending of the number of communities
  n_columns <- 3
  if (length(unique(tf_idf_table[Size_com >= threshold_com]$Com_wrap)) > 9) {
    n_columns <- 4
  }
  
  if (length(unique(tf_idf_table[Size_com >= threshold_com]$Com_wrap)) > 12) {
    n_columns <- 5
  }
  
  if (length(unique(tf_idf_table[Size_com >= threshold_com]$Com_wrap)) > 15) {
    n_columns <- 6
  }
  
  # plotting the graph
  tf_idf_plot <- ggplot(tf_idf_table[Size_com >= threshold_com], aes(reorder_within(term, count, color), count, fill = color)) +
    geom_bar(stat = "identity", alpha = .8, show.legend = FALSE) +
    labs(
      title = "Highest tf-idf",
      x = "Words", y = "tf-idf"
    ) +
    facet_wrap(~Com_wrap, ncol = n_columns, scales = "free") +
    scale_x_reordered() +
    scale_fill_identity() +
    theme(strip.text = element_text(size = size_title_wrap)) +
    coord_flip()
  
  list_return <- list("plot" = tf_idf_plot, "list_words" = tf_idf_table)
  return(list_return)
}
