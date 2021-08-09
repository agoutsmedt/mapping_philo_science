Script for building the different types of networks
================
Aurélien Goutsmedt
/ Last compiled on 2021-08-09

-   [1 What is this script for?](#what-is-this-script-for)
-   [2 Loading packages, paths and
    data](#loading-packages-paths-and-data)
    -   [2.1 External scripts](#external-scripts)
    -   [2.2 Loading Data](#loading-data)
-   [3 Bibliographic Coupling network](#bibliographic-coupling-network)
    -   [3.1 Building the network](#building-the-network)

# 1 What is this script for?

In this script, we build different networks (cocitation, coupling and
semantic similarity) for different subperiods. We use the Leiden
algorithm to identify communities, and calculate spatial coordinates
using Force Atlas 2 algorithm.

> WARNING: This script represents a first step of the project, and some
> processes have been improved (notably by the creation of new
> functions).

# 2 Loading packages, paths and data

## 2.1 External scripts

``` r
source("functions.R")
source("packages_and_paths.R")
```

## 2.2 Loading Data

``` r
direct_citations <- readRDS(paste0(data_path, "cleaned_direct_citations.rds")) %>% 
  as.data.table()
citing_articles <- readRDS(paste0(data_path, "cleaned_citing_articles.rds")) %>% 
  as.data.table()
```

# 3 Bibliographic Coupling network

The first step is to build the network from the citations data. Then we
will identify different communities through the Leiden algorithm, and
use the Force Atlas 2 algorithm to find spatial coordinates for each
citing document. Finally we will be able to project the graph.

## 3.1 Building the network

``` r
edges_threshold <- 2
nodes_coupling_threshold <- 0
```

We have to make first choices here. We have to decide if we take all the
citing articles, even if they are not cited at all by other documents in
our corpus. Here, we have decided to include all citing documents. For
the links between articles, we use a simple threshold method: we
consider a link between two articles important only if they share 2
references in common.

We have all the authors listed in the `citing_articles` table, so we
keep only the first authors for the coupling. We also keep only the
needed information.

``` r
nodes_coupling <- citing_articles[Seq_No == 1, c("OST_BK","new_ID_Ref","Annee_Bibliographique","Full_Name","revue","titre","Abstract")]
```

``` r
nb_cit <- direct_citations[, nb_cit := .N, by = "new_ID_Ref"]
nodes_coupling <- merge(nodes_coupling, unique(nb_cit[, c("new_ID_Ref", "nb_cit")]), by = "new_ID_Ref", all.x = "TRUE")

rm(nb_cit) # not useful anymore
# replacing NA value by 0 in nb_cit (it means that the nodes will be used in Leiden and Force Atlas, but
# not displayed in the graph)
nodes_coupling[is.na(nb_cit), ]$nb_cit <- 0
```

We will look into the distribution of citations in our corpus and will
reduce the number of nodes depending on the threshold
(`nodes_coupling_threshold`) we have chosen.

``` r
ggplot(nodes_coupling, aes(x = nb_cit)) +
  geom_boxplot(fill = "blue")
```

![](2_build_network_files/figure-gfm/r%20distrib_citations-1.png)<!-- -->

``` r
# reducing the number of nodes depending of the number of citations
nodes_coupling <- nodes_coupling[nb_cit >= nodes_coupling_threshold]
```

The next step just aim at securing that we only take citing documents in
the `direct_citations` that are also within the nodes. We also remove
articles without any references

``` r
edges_coupling <- unique(direct_citations[OST_BK %in% nodes_coupling$OST_BK])
if(length(edges_coupling$OST_BK) == length(direct_citations$OST_BK)){
  message("Everything was ok in the data extraction from the OST db: all the citing documents in the 
          direct citations table are also in the citing documents table.")
}
nodes_coupling <- nodes_coupling[OST_BK %in% edges_coupling$OST_BK]
```

We can now create the list of edges. We create two objects, with two
different methods (see the
[biblionetwork](https://agoutsmedt.github.io/biblionetwork/) package
documentation for more information):

-   the coupling angle value
-   the coupling strength value

``` r
edges_coupling_angle <- biblio_coupling(edges_coupling, source = "OST_BK", ref = "new_ID_Ref", weight_threshold = edges_threshold)
edges_coupling_strength <- coupling_strength(edges_coupling, source = "OST_BK", ref = "new_ID_Ref", weight_threshold = edges_threshold)
```

We remove the nodes that are not any more linked because they did not
pass the threshold with any other articles. We can now use the tidygraph
and the networkflow packages and build the two networks.

``` r
nodes_coupling_angle <- nodes_coupling[OST_BK %in% c(edges_coupling_angle$from,edges_coupling_angle$to)]
nodes_coupling_strength <- nodes_coupling[OST_BK %in% c(edges_coupling_strength$from,edges_coupling_strength$to)]

nodes_coupling_angle[, new_ID_Ref := as.character(OST_BK)]
nodes_coupling_strength[, new_ID_Ref := as.character(OST_BK)]

graph_coupling_angle <- tbl_main_component(edges_coupling_angle, nodes_coupling_angle, directed = FALSE, node_key = "OST_BK")
graph_coupling_strength <- tbl_main_component(edges_coupling_strength, nodes_coupling_strength, directed = FALSE, node_key = "OST_BK")
```

With the `edges_threshold` fixed at 2, we have removed 415 documents.
Before to move further, one small things we can do is to compare the
distribution of edges weights for the two methods.

``` r
compare_distrib <- data.table(angle_weight = edges_coupling_angle$weight,
                              strength_weight = edges_coupling_strength$weight)
compare_distrib <- pivot_longer(compare_distrib, cols = 1:2, names_to = "Measure", values_to = "Weight")

ggplot(compare_distrib, aes(Weight, after_stat(count), fill = Measure)) +
  geom_density(alpha = 0.3) +
  xlim(0, 0.35) +
  ggplot2::annotate("text",
                    x = 0.2,
                    y = 30000000,
                    label = paste0("max weight for coupling angle: ", max(edges_coupling_angle$weight))) +
  ggplot2::annotate("text",
                    x = 0.2,
                    y = 28000000,
                    label = paste0("max weight for coupling strength: ", max(edges_coupling_strength$weight)))
```

![](2_build_network_files/figure-gfm/r%20test_edges_distribution-1.png)<!-- -->

Regarding the distribution of the coupling strength, which is much less
widespread, it seems preferable, in a first time to opt for the coupling
angle measure, which will allow more difference in the weights of edges.

``` r
#Layout
#graph_before <- tbl %>% activate(nodes) %>% mutate(id=as.character(Id))
#nodes_before <- tbl %>% activate(nodes) %>% as.data.table()

if("x" %in% colnames(nodes_before)){
  graph_before <- graph_before %>% activate(nodes) %>% mutate(x = ifelse(is.na(x)==TRUE,sample(100:500, 1),x))
  graph_before <- graph_before %>% activate(nodes) %>% mutate(y = ifelse(is.na(y)==TRUE,sample(100:500, 1),y))
  graph_before <- graph_before %>% activate(nodes) %>% select(c(id,ID_Art,x,y,size))
} else{
  graph_before <- graph_before %>% activate(nodes) %>% select(c(id,ID_Art,size))
}

write.graph(graph = graph_before, file = 'tidy.graphml', format = 'graphml')
system("java -jar /projects/data/alexandre/GephiLayouts-1.0.jar forceatlas2 -i ./tidy.graphml -o  ./forceatlas2.graphml -threads 16 -maxiters 30000 -barneshut true -adjustsizes true -gravity 1")
gml <- read.graph("forceatlas2.graphml", format="graphml")
graph_after <- as_tbl_graph(gml)
graph_after <- graph_after %>% activate(nodes) %>% as.data.table() %>% .[,.(x,y,ID_Art)]

if("x" %in% colnames(nodes_before)){
  tbl <- tbl %>% activate(nodes) %>% select(-c(x,y)) %>% left_join(graph_after)
} else{
  tbl <- tbl %>% activate(nodes) %>% left_join(graph_after)
}
```
