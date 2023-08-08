library(tidyverse)

bing_coho_up <- read_rds("data/bing_coho_up.rds") |>
  mutate(sex = ifelse(type == "female", "female", "male"),
         age = ifelse(type == "jack", 2, 3)) 

bing_age <- bing_coho_up |>
  group_by(year, age) |>
  summarize(count = sum(count),
            .groups = "drop_last") |>
  mutate(pct = count / sum(count)) |>
  ungroup()



bing_age |> filter(age == 2) |>
ggplot(aes(x = as.integer(year), y = pct)) +
  geom_line() +
  geom_hline(aes(yintercept = mean(pct)), color = "blue") +
  geom_hline(aes(yintercept = weighted.mean(pct, count)), color = "red")

bing_age |> 
  group_by(year) |>
  summarize(count = sum(count),
            pct = min(pct)) |>
  ggplot(aes(x = count, y = pct)) +
  geom_point()

bing_age |>
  select(year, age, pct) |>
  pivot_wider(id_cols = year,
              names_from = age,
              names_prefix = "age_",
              values_from = pct) |>
  mutate(log_ratio = log(age_3) / log(age_2)) |>
  summarize(age_2 = mean(age_2),
            age_3 = mean(age_3),
            sd_l2 = sd(log(age_2)),
            sd_alr = sd(log_ratio))
