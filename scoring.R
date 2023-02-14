

scores <- neon4cast::combined_scores(theme="aquatics", collect=FALSE)

model_ids <- c("xgboost_temp_oxygen_parallel", "xgboost_temp_oxygen_sequential", "climatology", "flareGLM")

df <- scores %>% 
  dplyr::filter(model_id %in% model_ids) %>% 
  dplyr::filter(reference_datetime > "2023-02-01") %>% 
  dplyr::collect()
