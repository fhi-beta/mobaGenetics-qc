library(janitor)
library(glue)
library(tidyr)
library(dplyr)
library(stringr)
library(knitr)


debug <- F

if (debug) {
  
  args <- c("/mnt/archive3/phasing_test/expected_all_relations", "/mnt/archive3/phasing_test/new_imputation_batches", "/mnt/archive3/phasing_test/batch_shuffle.md",  "/mnt/archive3/phasing_test/batch_movements", "/mnt/archive3/phasing_test/batch_movements_trios"))
  
} else {
  args <- commandArgs(TRUE)
}


rel_file <- args[1]
new_batches_trunk <- args[2]
md_file <- args[3]
movements_file <- args[4]
trios_file <- args[5]

rel <- read.table(rel_file, header = T)
rel$orig_batch <- rel$iid_batch
rel$orig_parents_in_batch <- rel$parents_in_batch


problem_children <- subset(rel, (!is.na(pat) | !is.na(mat)) & parents_in_batch == 0)
to_move <- subset(problem_children, shared_chips > 0)

to_move$move_from <- to_move$iid_batch
to_move$move_to <- NA

to_move <- to_move %>% mutate(move_to = ifelse(!is.na(mat_chip) & iid_chip == mat_chip,
                                   mat_batch, 
                                   pat_batch))

to_move$moved <- 1


updated_rel <- rel %>% left_join(to_move %>% select(iid, move_from, move_to, moved), by = "iid") %>% mutate(iid_batch = ifelse(!is.na(move_to), move_to, iid_batch))

updated_rel <- updated_rel %>%
  mutate(parents_in_batch = ifelse(
    !is.na(mat_batch) & !is.na(pat_batch) & iid_batch == mat_batch & iid_batch == pat_batch, 2,
    ifelse(
      (!is.na(mat_batch) & iid_batch == mat_batch) | (!is.na(pat_batch) & iid_batch == pat_batch), 1, 
      0
    )
  ))



batches <- unique(updated_rel$iid_batch)
updated_trios <- subset(updated_rel, !is.na(pat) & !is.na(mat))

for (batch in batches) {
  samples_file <- paste0(new_batches_trunk, ".imputation.", batch)
  samples <- subset(updated_rel, imputation_batch == batch) %>% select(iid)
  samples <- unique(samples)
  write.table(samples, samples_file, row.names = F, col.names = F, quote=F)
}

write.table(x = updated_rel, file = movements_file, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(x = updated_trios, file = trios_file, col.names = F, row.names = F, quote = F, sep = "\t")



# problem_children_updated <- subset(updated_rel, parents>0 & shared_chips == 0)
# no_problem_children_updated <-subset(updated_rel, parents>0 & shared_chips > 0)
# problem_fathers <- unique(na.omit(problem_children_updated$pat))
# problem_fathers_without_unproblematic_children <- problem_fathers[!(problem_fathers %in% no_problem_children_updated$pat)]
# problem_mothers <- unique(na.omit(problem_children_updated$mat))
# problem_mothers_without_unproblematic_children <- problem_mothers[!(problem_mothers %in% no_problem_children_updated$mat)]

# updated_rel <- updated_rel %>% mutate(iid_batch = ifelse(iid %in% problem_children_updated$iid, "problem", iid_batch))
# updated_rel <- updated_rel %>% mutate(iid_batch = ifelse(iid %in% problem_fathers_without_unproblematic_children | iid %in% problem_mothers_without_unproblematic_children, "problem", iid_batch))

# updated_rel <- updated_rel %>% left_join(imputation_merge %>% select (iid_batch = orig_batch, imputation_batch = new_batch), by = "iid_batch")


# # Move problem to separate batch (with problem parents) for phasing
# for (batch in batches) {
#   samples_file <- paste0(new_batches_trunk, ".phasing.", batch)
#   samples <- subset(updated_rel, iid_batch == batch) %>% select(iid)
#   samples <- unique(samples)
#   write.table(samples, samples_file, row.names = F, col.names = F, quote=F)
# }




# problem_batch_file <- paste0(new_batches_trunk, ".phasing.problem")

# problem_children_char <- as.character(problem_children_updated$iid)
# problem_fathers_char <- as.character(problem_fathers)
# problem_mothers_char <- as.character(problem_mothers)
# problem_fathers_without_unproblematic_children_char <- as.character(problem_fathers_without_unproblematic_children)
# problem_mothers_without_unproblematic_children_char <- as.character(problem_mothers_without_unproblematic_children)


# writeLines(problem_children_char, problem_batch_file)
# writeLines(problem_children_char, problem_children_file)

# con <- file(problem_batch_file, open = "a")
# writeLines(problem_fathers_char, con)
# close(con)

# con <- file(problem_batch_file, open = "a")
# writeLines(problem_mothers_char, con)
# close(con)

# con <- file(problem_children_file, open = "a")
# writeLines(problem_fathers_without_unproblematic_children_char, con)
# close(con)

# con <- file(problem_children_file, open = "a")
# writeLines(problem_mothers_without_unproblematic_children_char, con)
# close(con)






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

