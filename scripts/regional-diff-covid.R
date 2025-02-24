library(tidyverse)
library(cowplot)
library(scoringutils)

covid <- read_csv("data/weekly-covid-us.csv")

res <- read_csv("~/Downloads/coviddatelocation.csv") |> 
    select(model:date, r_wis)

us <- load_us_graph(flu) |> # for HHS region codes
    st_drop_geometry() |> 
    mutate(region=ifelse(region == 9, 3, region)) |> 
    mutate(region=factor(region, labels=c("Northeast", "Midwest", "South", "West/Pacific")))


reg_diff <- res |> 
    filter(model != "baseline_covid") |> 
    left_join(us, by=join_by(location == state)) |> 
    group_by(date, region) |> 
    summarize(reg_rwis=mean(r_wis), .groups="drop") |> 
    pivot_wider(names_from=region, values_from=reg_rwis) |> 
    mutate(across(-date, ~. - South)) |> 
    pivot_longer(-date) |> 
    filter(name != "South")

#66C5CC,#F6CF71,#F89C74,#DCB0F2,#87C55F,#9EB9F3,#FE88B1,#C9DB74,#8BE0A4,#B497E7,#D3B484,#B3B3B3

p1 <- ggplot(reg_diff, aes(date, value, col=name)) +
    geom_hline(yintercept=0, col="gray70") +
    geom_line(linewidth=1.01) +
    scale_color_manual(values=c("#87C55F","#FE88B1","#B497E7")) +
    scale_x_date(date_breaks="3 months", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="\nrWIS difference", col=NULL) +
    theme_half_open() +
    theme(axis.text.x=element_text(size=rel(0.9)), legend.position="none")

p2 <- reg_diff |> 
    group_by(name) |> 
    mutate(cum_diff=cumsum(value)) |> 
    ggplot(aes(date, cum_diff, col=name)) +
    geom_line(linewidth=1.01) +
    scale_color_manual(values=c("#87C55F","#FE88B1","#B497E7")) +
    scale_x_date(date_breaks="3 months", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="\nCumulative difference", col=NULL) +
    theme_half_open() +
    theme(axis.text.x=element_text(size=rel(0.9)), legend.position="none")

covid_reg <- covid |> 
    left_join(us, by=join_by(location == state)) |> 
    filter(!is.na(region), date >= "2022-09-17") |> 
    group_by(date, region) |> 
    summarize(reg_rate=mean(weekly_rate))

covid_nat <- covid |> 
    filter(location == "US", date >= "2022-09-17") |> 
    slice(seq(2, n(), 2)) |> 
    mutate(location="National")

p3 <- ggplot(covid_reg, aes(date, reg_rate)) +
    geom_line(aes(col=fct_relevel(region, "Midwest", "Northeast", "West/Pacific")), linewidth=1.01) +
    geom_point(
        aes(y=weekly_rate, col=location), data=covid_nat,
        alpha=0.8, shape=1
    ) +
    scale_y_continuous(transform="sqrt") +
    scale_color_manual(values=c("#87C55F", "#FE88B1", "#B497E7", "#66C5CC", "gray30")) +
    scale_x_date(date_breaks="3 months", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="\nAdmits per 100k", col=NULL) +
    theme_half_open() +
    theme(axis.text.x=element_text(size=rel(0.9)))

plot_grid(
    plot_grid(p3 + theme(legend.position="none"), p1, p2, nrow=3, labels="AUTO", align="v", axis="l"),
    get_legend(p3),
    nrow=1, rel_widths=c(1, 0.3)
)

# plot_grid(
#     p3,
#     plot_grid(
#         plot_grid(p1, p2, nrow=2, labels=c("B", "C"), align="v", axis="l"),
#         get_legend(p1 + theme(legend.position="right")),
#         rel_widths=c(1, 0.3)   
#     ),
#     labels=c("A", ""), nrow=2, rel_heights=c(0.36, 2/3)
# )

ggsave("figs/regional-diff-covid4.pdf", width=7, height=6.3)
