biogrid_data <- read.delim("/home/francescoc/Downloads/BIOGRID-ALL-4.4.241.tab3.txt", header = TRUE, stringsAsFactors = FALSE)

biogrid_data <- biogrid_data[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B", "Score","Experimental.System.Type")]
names(biogrid_data) <- c("Interactor_A", "Interactor_B", "Score", "type")

biogrid_data <- biogrid_data %>%
  mutate(Score = as.numeric(Score)) %>%
  filter(type == "physical") %>%       
  filter(!is.na(Score)) %>%            
  filter(abs(Score) > 300)

write.table(biogrid_data, "/home/francescoc/Desktop/GRN_project/data/biogrid_physical_s400.txt", sep = "\t", quote = F, col.names = T, row.names = T)
