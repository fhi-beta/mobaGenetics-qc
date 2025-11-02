library(janitor)
library(glue)
library(tidyr)
library(dplyr)
library(stringr)
library(knitr)


debug <- F

if (debug) {
  
  args <- c("/mnt/archive3/phasing_test/all_relations", "/mnt/archive3/phasing_test/new_batches", "/mnt/archive3/phasing_test/batch_shuffle.md")
  
} else {
  args <- commandArgs(TRUE)
}


rel_file <- args[1]
new_batches_trunk <- args[2]
md_file <- args[3]

rel <- read.table(rel_file, header = T)
rel$orig_batch <- rel$iid_batch
rel$orig_parents_in_batch <- rel$parents_in_batch



problem_children <- subset(rel, (pat != "0" | mat != "0") & parents_in_batch == 0)
problem_children$move_from <- problem_children$iid_batch
problem_children$move_to <- ifelse(!is.na(problem_children$mat_batch), 
                                   problem_children$mat_batch, 
                                   problem_children$pat_batch)
problem_children$moved <- 1

updated_rel <- rel %>% left_join(problem_children %>% select(iid, move_from, move_to, moved), by = "iid") %>% mutate(iid_batch = ifelse(!is.na(move_to), move_to, iid_batch))

updated_rel <- updated_rel %>%
  mutate(parents_in_batch = ifelse(
    !is.na(mat_batch) & !is.na(pat_batch) & iid_batch == mat_batch & iid_batch == pat_batch, 2,
    ifelse(
      (!is.na(mat_batch) & iid_batch == mat_batch) | (!is.na(pat_batch) & iid_batch == pat_batch), 1, 
      0
    )
  ))


batches <- unique(updated_rel$iid_batch)
for (batch in batches) {
  samples_file <- paste0(new_batches_trunk, ".batch_", batch)
  samples <- subset(updated_rel, iid_batch == batch) %>% select(iid)
  write.table(samples, samples_file, row.names = F, col.names = F, quote=F)
}



# Calculate Movement Matrix
movement_matrix <- subset(updated_rel, moved == 1) %>%
  select(move_from, move_to) %>%
  group_by(move_from, move_to) %>%
  summarize(count = n(), .groups = 'drop') %>%
  spread(key = move_to, value = count, fill = 0)

# Calculate Batch Summary
orig_count <- updated_rel %>%
  group_by(orig_batch) %>%
  summarize(original_count = n(), .groups = 'drop')

movements_out <- updated_rel %>%
  group_by(move_from) %>%
  summarize(out_count = n(), .groups = 'drop')

movements_in <- updated_rel %>%
  group_by(move_to) %>%
  summarize(in_count = n(), .groups = 'drop')

batch_summary <- orig_count %>%
  left_join(movements_out, by = c("orig_batch" = "move_from")) %>%
  left_join(movements_in, by = c("orig_batch" = "move_to")) %>%
  mutate(out_count = replace_na(out_count, 0), in_count = replace_na(in_count, 0)) %>%
  mutate(new_count = original_count - out_count + in_count)

# Format tables as markdown and add captions
batch_summary_md <- kable(batch_summary, format = "markdown", caption = "Batch Summary", 
                          col.names = c("Batch", "Original Count", "Movements Out", "Movements In", "New Count"))

movement_matrix_md <- kable(movement_matrix, format = "markdown", caption = "Movement Matrix")

write(
  x = batch_summary_md,
  file = md_file,
  append = F
)

write(
  x = movement_matrix_md,
  file = md_file,
  append = T
)

