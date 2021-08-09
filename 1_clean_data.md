Script for cleaning the data
================
Aurélien Goutsmedt
/ Last compiled on 2021-08-09

-   [1 What is this script for?](#what-is-this-script-for)
-   [2 Loading packages, paths and
    data](#loading-packages-paths-and-data)
    -   [2.1 External scripts](#external-scripts)
    -   [2.2 Loading Data](#loading-data)
-   [3 Cleaning the data](#cleaning-the-data)
    -   [3.1 Cleaning titles and authors](#cleaning-titles-and-authors)
    -   [3.2 Building the identifiers](#building-the-identifiers)
    -   [3.3 Cleaning the nodes with new
        ID](#cleaning-the-nodes-with-new-id)

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
direct_citations <- read_excel(paste0(data_path, "Philo_Sciences_Couplage-Cocitation_2010-2020.xlsx")) %>% 
  as.data.table()
citing_articles <- read_excel(paste0(data_path, "Philo_des_Sciences_2010-2020_Gender_VL.xlsx")) %>% 
  as.data.table()
```

# 3 Cleaning the data

In the OST database, only the articles, books and book chapters that are
listed in the citing documents table have a unique identifier in the
cited documents table. It means that many books have no identifier (more
precisely, there is one identifier per citation, meaning that the same
book will have as many identifiers as the number of times it is cited).
Thus, they will be excluded from the coupling and co-citation analysis,
which corresponds to a large loss of information (40% of citations have
no proper identifier).

This part of the script aims at building consistent identifiers for
those references which lacks one.

## 3.1 Cleaning titles and authors

The first step is to normalized titles and authors to allow the creation
of ID depending on titles and authors.

``` r
direct_citations[, Cited_Work := toupper(Cited_Work)]

# small cleaning for problems with punctuations
direct_citations[str_which(Cited_Work, "[A-z]- [A-z]|[A-z] -[A-z]"), Cited_Work := str_replace_all(Cited_Work, " -", "-")]

# cleaning authors
direct_citations[, Cited_Author := toupper(Cited_Author)]
direct_citations[, Cited_Author := str_replace_all(Cited_Author, ", ","-")]
direct_citations[, Cited_Author := str_replace_all(Cited_Author, "\\. ","")] # Removing point after Initials
direct_citations[, Cited_Author := str_replace(Cited_Author, "\\.","")] # Removing the last point (which has no space after it)
```

We create a new column to save the original column, just in a case of
problems when we clean the names.

``` r
direct_citations[, cleaned_cited_author := str_replace_all(Cited_Author, " |\\.","-")]
direct_citations[, cleaned_cited_author := str_replace_all(cleaned_cited_author, ",","")]
direct_citations[, cleaned_cited_author := str_replace(cleaned_cited_author, "-$","")] # Removing the last hyphen (which has no space after it)
direct_citations[, cleaned_cited_author := str_replace(cleaned_cited_author, "--","-")]

# We now need to transform in initials all the first names that remains after the dash.

direct_citations[, post_dash := str_extract(cleaned_cited_author, "\\-.*")]
direct_citations[, n_dash := str_count(post_dash, "-")] # will serve for distinguish between the different cases

direct_citations[n_dash == 1 &
                   str_detect(post_dash, "^-") &
                   str_length(post_dash) > 4 &
                   ! str_detect(post_dash, "GOVERNMENT|INTERNATIONAL|\\(|FOUNDATION"),
                 post_dash := str_sub(post_dash, 1, 2)]
```

When we have two dashes, we want to keep the first letter of the first
and the second name.

``` r
direct_citations[n_dash == 2, `:=` (initial_1 = str_sub(str_remove(str_extract(post_dash, "-[A-Z]+-"), "-$"), 1, 2),
                                    initial_2 = str_sub(str_remove(str_extract(post_dash, "-[A-Z]+$"), "^-"), 1, 1))]
direct_citations[n_dash == 2, post_dash := paste0(initial_1,initial_2)]
```

We now transformed `cleaned_cited_author` using the `post_dash` column
that we have cleaned.

``` r
direct_citations[, cleaned_cited_author := paste0(str_remove(cleaned_cited_author, "\\-.*"),post_dash)]
```

Some special cases remains (around 3000, gathering things with one or
two dash that we have no cleaned properly, or things with three dashes).
But as it is minor cases (some names of institutions, some unusual
names, etc.) it should not impact what this cleaning aimed at:
standardizing names to build common identifiers.

## 3.2 Building the identifiers

The strategy is to consider that:

1.  ref with the same title, same date and the same author before the
    dash are the same ref.
2.  ref with same author and same title are the same ref.
3.  If same author, same date, and similarities in the title, same ref.

``` r
# Copying id
direct_citations[, cleaned_ID_Ref := UID_Ref]

# That is the dt we will use for each method
to_identify <- direct_citations[str_detect(UID_Ref, "\\.[0-9]{1,3}$"), c("UID_Ref","cleaned_cited_author","Cited_Work","Year")] %>% 
  .[order(UID_Ref)] %>% 
  unique()
```

The goal is to build a dt with former and new ID after having
implemented the three different methods. \#\#\# Method 1: ref with the
same title, same date and the same author before the dash

``` r
method_1 <- to_identify[, just_surname := str_trim(str_remove(str_extract(cleaned_cited_author, "[A-Z]*-"), "-$"))] %>% 
  .[! is.na(just_surname) & ! just_surname == "", c("UID_Ref","just_surname","Cited_Work","Year")] 
doublons <- which(duplicated(method_1[, c("just_surname","Cited_Work","Year")]))
new_ID <- method_1[-doublons]
setnames(new_ID, "UID_Ref","new_ID")
method_1 <- merge(method_1, new_ID, by = c("just_surname","Cited_Work","Year"))
method_1 <- method_1 %>% 
  select(UID_Ref, new_ID) %>% 
  filter(UID_Ref != new_ID) %>% 
  unique()
```

We build a simple test to verify that a former ID is attributed to only
one new ID.

``` r
test_dt <- method_1 %>% 
  group_by(UID_Ref) %>% 
  mutate(test = n()) %>% 
  arrange(UID_Ref) %>% 
  filter(test > 1)

if(length(test_dt$new_ID != 0)){
  warning("There is a problem here, with different new_ID for a unique UID_Ref")
  doublons <- which(duplicated(method_1$UID_Ref))
  method_1 <- method_1[-doublons]
} else {
  message("Everything ok: a unique new_ID for each UID_Ref")
}
```

Method 2: ref with same author and same title

``` r
method_2 <- to_identify[, c("UID_Ref","cleaned_cited_author","Cited_Work")] %>% 
  unique()
doublons <- which(duplicated(method_2[, c("cleaned_cited_author","Cited_Work")]))
new_ID <- method_2[-doublons]
setnames(new_ID, "UID_Ref","new_ID")
method_2 <- merge(method_2, new_ID, by = c("cleaned_cited_author","Cited_Work"))
method_2 <- method_2 %>% 
  select(UID_Ref, new_ID) %>% 
  filter(UID_Ref != new_ID) %>% 
  arrange(UID_Ref) %>% 
  unique()

test_dt <- method_2 %>% 
  group_by(UID_Ref) %>% 
  mutate(test = n()) %>% 
  filter(test > 1)

if(length(test_dt$new_ID != 0)){
  warning("There is a problem here, with different new_ID for a unique UID_Ref")
  doublons <- which(duplicated(method_2$UID_Ref))
  method_2 <- method_2[-doublons]
} else {
  message("Everything ok: a unique new_ID for each UID_Ref")
}
```

Method 3: same author, same date, and similarities in the title

The strategy here is to isolates the words of a title, and to remove
stop words before reconstituting titles. We keep only the titles with at
least three words because it would be too risky otherwise.

``` r
reduce_title <- to_identify %>% 
  unnest_tokens(word, Cited_Work) %>% 
  anti_join(stop_words) %>% 
  group_by(UID_Ref) %>% 
  mutate(n_word = n()) %>% 
  filter(n_word > 2)

reduce_title <- reduce_title %>% 
  group_by(UID_Ref) %>% 
  summarize(reduce_title = paste(word, collapse = " ")) %>% 
  ungroup()

method_3 <- merge(to_identify, reduce_title, by = "UID_Ref")
method_3 <- method_3[order(UID_Ref)]

doublons <- which(duplicated(method_3[, c("cleaned_cited_author","Year", "reduce_title")]))
new_ID <- method_3[-doublons]
setnames(new_ID, "UID_Ref","new_ID")
method_3 <- merge(method_3, new_ID, by = c("cleaned_cited_author","Year", "reduce_title"))
method_3 <- method_3 %>% 
  select(UID_Ref, new_ID) %>% 
  filter(UID_Ref != new_ID) %>%
  arrange(UID_Ref) %>% 
  unique()

test_dt <- method_3 %>% 
  group_by(UID_Ref) %>% 
  mutate(test = n()) %>% 
  filter(test > 1)

if(length(test_dt$new_ID != 0)){
  warning("There is a problem here, with different new_ID for a unique UID_Ref")
  doublons <- which(duplicated(method_3$UID_Ref))
  method_3 <- method_3[-doublons]
} else {
  message("Everything ok: a unique new_ID for each UID_Ref")
}
```

We can now bind the results for the three methods, and check that there
is no inconsistency, like 1) different `new_id` attributed to the same
`UID_Ref` or 2) a `new_id` that is also a `UID_Ref` in the table.

``` r
bind_method <- rbind(method_1,method_2,method_3) %>% 
  unique() %>% 
  arrange(UID_Ref, new_ID)
```

The strategy for problem 1) is simply to remove the doublons in
`UID_Ref`, as we assume that any `new_ID` is ok but we just have to keep
one.

``` r
doublons <- which(duplicated(bind_method$UID_Ref))
bind_method <- bind_method[-doublons]
```

The strategy is to identify the `new_ID` which are in `UID_Ref` and to
replace this `new_ID` by the `new_ID` associated to the `UID_Ref`. This
strategy have to be replicated until there is no more overlapping.

``` r
for(i in 1:100) {
  problem_2 <- bind_method[UID_Ref %in% bind_method$new_ID]
  
  if (length(problem_2$UID_Ref) == 0) {
    break
  } else {
    warning(paste0(
      "round ",
      i,
      ": there are ",
      length(problem_2$UID_Ref),
      " ID that are both in UID_Ref and new_ID"
    ))
    
    setnames(problem_2, "new_ID", "new_new_ID")
    bind_method <- merge(
      bind_method,
      problem_2,
      by.y = "UID_Ref",
      by.x = "new_ID",
      all.x = TRUE
    )
    bind_method[!is.na(new_new_ID)]$new_ID <-
      bind_method[!is.na(new_new_ID)]$new_new_ID
    bind_method <- bind_method[, c("UID_Ref", "new_ID")] %>%
      .[order(UID_Ref, new_ID)] %>%
      unique()
  }
  
}
```

We can now cleaned ID in the main file, and rebuild a totally new ID
from that to simplify merging later (current ID are in character format
because of characters in the ID, and it would improve effiency to have
ID in an integer format).

``` r
direct_citations <- merge(direct_citations, bind_method, by = "UID_Ref", all.x = TRUE)
direct_citations[! is.na (new_ID)]$cleaned_ID_Ref <- direct_citations[! is.na (new_ID)]$new_ID 
new_ID <- data.table(cleaned_ID_Ref = unique(direct_citations$cleaned_ID_Ref), 
                     new_ID_Ref = 1:length(unique(direct_citations$cleaned_ID_Ref)))
direct_citations <- merge(direct_citations, new_ID, by = "cleaned_ID_Ref")
direct_citations <- direct_citations[, -"new_ID"]
```

## 3.3 Cleaning the nodes with new ID

We can now integrate the new IDs to the `citing_articles` file

``` r
new_ID <- direct_citations[, c("UID_Ref", "new_ID_Ref")] %>% 
  unique()
citing_articles <- merge(citing_articles, new_ID, by.x = "UID", by.y = "UID_Ref", all.x = TRUE)
```

We can now save the cleaned data to be load in following scripts

``` r
saveRDS(direct_citations, paste0(data_path, "cleaned_direct_citations.rds"))
saveRDS(citing_articles, paste0(data_path, "cleaned_citing_articles.rds"))
```
