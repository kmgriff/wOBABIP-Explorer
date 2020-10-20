##########################################################
# bipexplorer.R
# Imports MLB Statcast data and precomputed regression objects from local storage, slices
# them into useful subsets, then defines an RShiny server to create interactive plots.
#
# The full data plot generates wOBABIP heatmaps over all of the data between 2017-2020.
# Interactive sliders can be used to limit the data observed.  Tip: Try setting it to
# "Launch speed", and limit the plot to 0-50 MPH.  Do players tend to bunt toward first base,
# or third base?
#
# The Team Explorer shows the plot for specific teams, and can be further subsetted by
# year.  If they are the hitting team, wOBABIP shows how they perform, while in fielding mode
# it shows every other team's performance against them.  Checking "Difference plot" will show
# how their offense / defense compares against the total dataset.  Checking "Away games only"
# reduces the impact of different stadium dimensions.
#
# TODO: Checking "Difference" will show how their offense / defense compares against the
# average of all the teams' regressions for the same subsets.
##########################################################

library(tidyverse)
library(lubridate)
library(np)
library(sp)
library(metR)
library(geometry)
library(shiny)
library(shinythemes)
library(RColorBrewer)
library(plotly)

totaldat <- readRDS("data/all_mlb.RDS")
totalwobareg <- readRDS("data/gridwobafull.RDS")
totalwobareg2d <-
  lapply(as.character(2017:2020), function(x) {
    readRDS(sprintf("data/gridwobafull2d_%s.RDS", x))
  })
names(totalwobareg2d) <- 2017:2020
totaldensreg <- readRDS("data/griddensfull.RDS")
totaldensreg2d <-
  lapply(as.character(2017:2020), function(x) {
    readRDS(sprintf("data/griddensfull2d_%s.RDS", x))
  })
names(totaldensreg2d) <- 2017:2020
totalgrid <- readRDS("data/gridfull.RDS")
totalgrid2d <-
  lapply(as.character(2017:2020), function(x) {
    readRDS(sprintf("data/gridfull2d_%s.RDS", x))
  })
names(totalgrid2d) <- 2017:2020

teamids <- totaldat$home_team %>% unique() %>% sort()



totalwobadens <- totalgrid %>%
  mutate(wOBABIP = totalwobareg$mean, density = totaldensreg$dens)
totalwobadens2d <- map(as.character(2017:2020), ~{
  data.frame(get(.x, totalgrid2d)) %>% mutate(
    wOBABIP = get(.x, totalwobareg2d)$mean,
    density = get(.x, totaldensreg2d)$dens
    )
  })
names(totalwobadens2d) <- 2017:2020



calcTeamWoba <- function (team, awayOnly, year, fieldPosition, difference) {
  away <- ifelse(awayOnly, "away", "all")
  
  if(!difference) {
    teamfile <- readRDS(str_glue("data/reg2d_{team}_{fieldPosition}_{away}_{year}.RDS"))
    teamfile$reg <- teamfile$reg %>% rename(wOBABIP = mean, density = dens)
    list(teamfile$reg, teamfile$n)
    
  } else {
    allfile <- map(teamids, ~readRDS(str_glue("data/reg2d_{.x}_{fieldPosition}_{away}_{year}.RDS")))
    teamfile <- allfile[[match(team, teamids)]]
    
    regs <- bind_rows(map(allfile, ~.x$reg %>% mutate(n = .x$n))) %>%
      group_by(launch_angle, launch_speed) %>%
      summarize(mean = weighted.mean(mean, dens),
                dens = weighted.mean(dens, n)) %>%
      ungroup()
    
    diffreg <- inner_join(teamfile$reg, regs, by = c("launch_angle", "launch_speed")) %>%
      mutate(wOBABIP = mean.x - mean.y,
             density = dens.x - dens.y)
    
    list(diffreg, teamfile$n)
  }
  # 
  # subset <- totaldat %>%
  #   filter(year(game_date) == year)
  # 
  # if (fieldPosition == "hitting") {
  #   subset <- subset %>% filter(hitting_team == team)
  # } else {
  #   subset <- subset %>% filter(fielding_team == team)
  # }
  # 
  # if (awayOnly) {
  #   subset <- subset %>% filter(away_team == team)
  # }
  # 
  # txdat <-
  #   data.frame(x = subset$launch_angle,
  #              y = subset$launch_speed)
  # 
  # totalgrid <- get(as.character(year), totalgrid2d)
  # 
  # hull <- convhulln(as.matrix(txdat))
  # grid_inhull <-
  #   totalgrid[inhulln(hull, as.matrix(totalgrid)), ]
  # 
  # reg <- npreg(
  #   bws = c(5, 40),
  #   bwtype = "adaptive_nn",
  #   exdat = grid_inhull,
  #   txdat = txdat,
  #   tydat = as.numeric(subset$woba_value)
  # )
  # 
  # density <- grid_inhull %>%
  #   map(~ {
  #     attributes(.x) <- NULL
  #     .x
  #   }) %>%
  #   as.data.frame() %>%
  #   npudens(
  #     bws = c(5, 40),
  #     bwtype = "adaptive_nn",
  #     edat = .,
  #     tdat = txdat
  #   )
  # 
  # out <-
  #   grid_inhull %>% mutate(wOBABIP = reg$mean, density = density$dens)
  # 
  # return(list(out, length(subset[[1]])))
}


plotWobaMargin <- function (l) {
  dat <- l[[1]]
  n <- l[[2]]
  
  breaks <- seq(0, 2, length.out = 9)
  dims <- names(dat)[1:2] %>%
    map(~case_when(
      .x == "launch_angle" ~ "Vertical angle",
      .x == "launch_speed" ~ "Launch speed",
      .x == "spray_angle" ~ "Horizontal angle"
      ))
  
  p <- ggplot(dat %>% 
                rename_all(~c("x", "y", "z")),
              aes(x = x, y = y, z = z)) +
      geom_contour_filled(breaks = breaks) +
      geom_contour_tanaka(breaks = breaks) +
      scale_fill_manual(values = brewer.pal(9, "Greens"), drop = FALSE)
    
    if (!"Vertical angle" %in% dims) {
      p <- p + coord_polar(theta = "y",
                           start = pi / 2,
                           direction = -1) +
        scale_y_continuous(
          breaks = seq(45, 135, by = 10),
          minor_breaks = NULL,
          limits = c(-180, 180)
        ) +
        scale_x_continuous(breaks = seq(0, 120, by = 10),
                           limits = c(0, 120))
    }
    else if (!"Horizontal angle" %in% dims) {
      p <- p + coord_polar(theta = "y",
                           start = pi / 2,
                           direction = -1) +
        scale_y_continuous(
          breaks = seq(-90, 90, by = 10),
          minor_breaks = NULL,
          limits = c(-180, 180)
        ) +
        scale_x_continuous(breaks = seq(0, 120, by = 10),
                           limits = c(0, 120))
    }
    else {
      p <- p + coord_cartesian() +
        scale_y_continuous(
          breaks = seq(-90, 90, by = 10),
          minor_breaks = NULL,
          limits = c(-90, 90)
        ) +
        scale_x_continuous(breaks = seq(45, 135, by = 10),
                           limits = c(45, 135))
    }
    p + labs(
      x = dims[1],
      y = dims[2],
      fill = "wOBABIP",
      title = "wOBABIP",
      caption = paste0("N = ", n)
    )
}

plotWobaTeam <- function (l, difference) {
  dat <- l[[1]]
  n <- l[[2]]
  
  if (difference) {
    colorpal = brewer.pal(11, "PiYG")
    breaks = seq(-.5, .5, length.out = 11)
  } else {
    colorpal = brewer.pal(9, "Greens")
    breaks <- seq(0, 2, length.out = 9)
  }
  
  ggplot(dat, aes(x = launch_speed, y = launch_angle, z = wOBABIP)) +
      geom_contour_filled(breaks = breaks) +
      geom_contour_tanaka(breaks = breaks) +
      scale_fill_manual(values = colorpal, drop = FALSE) +
      coord_polar(theta = "y",
                  start = pi / 2,
                  direction = -1) +
      scale_y_continuous(
        breaks = seq(-90, 90, by = 10),
        minor_breaks = NULL,
        limits = c(-180, 180)
      ) +
      scale_x_continuous(breaks = seq(0, 120, by = 10),
                         limits = c(0, 120)) +
      labs(
        x = "Launch speed",
        y = "Launch angle",
        fill = "wOBABIP",
        title = "wOBABIP",
        caption = paste0("N = ", n)
      )
}

# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
                navbarPage(
                  "BIP WOBA Viewer",
                  tabPanel("Full dataset",
                           sidebarLayout(
                             sidebarPanel(
                               # Select dimension to flatten
                               selectInput(
                                 inputId = "dim",
                                 label = strong("Select dimension for windowed average"),
                                 choices = c(
                                   "Vertical angle" = "launch_angle",
                                   "Horizontal angle" = "spray_angle",
                                   "Launch speed" = "launch_speed"
                                 )
                               ),
                               
                               
                               # Select range for marginal average
                               sliderInput(
                                 inputId = "range",
                                 label = strong("Select range for windowed average"),
                                 min = -180,
                                 max = 180,
                                 value = c(90, 90)
                               ),
                             ),
                             
                             mainPanel(plotOutput(outputId = "wobaplotfull", height = "800px"))
                           ),),
                  tabPanel("Team Explorer",
                           sidebarLayout(
                             sidebarPanel(
                               selectInput(
                                 inputId = "year",
                                 label = "Year",
                                 choices = 2017:2020
                               ),
                               selectInput(
                                 inputId = "team",
                                 label = "Team",
                                 choices = teamids
                               ),
                               selectInput(
                                 inputId = "fieldposition",
                                 label = "Field position",
                                 choices = c("Hitting" = "hitting",
                                             "Fielding" = "fielding")
                               ),
                               checkboxInput(inputId = "awayonly",
                                             label = "Away games only"),
                               checkboxInput(inputId = "difference",
                                             label = "Difference plot")
                             ),
                             
                             mainPanel(plotOutput(outputId = "wobaplotteam", height = "800px"))
                           ))
                ))

server <- function(input, output, session) {
  dims <- reactive({
    req(input$dim)
    setdiff(c("launch_speed", "spray_angle", "launch_angle"),
            c(input$dim))
  })
  
  dimRange <- reactive({
    req(input$dim)
    range(totalwobadens[[input$dim]])
  })
  
  # Subset and average data
  wobamarginalized <- reactive({
    req(input$dim, input$range)
    n <- totaldat %>%
      filter_at(input$dim, ~.x >= input$range[1] & .x <= input$range[2]) %>%
      {length(.[[1]])}
    totalwobadens %>%
      filter_at(input$dim, ~.x >= input$range[1] & .x <= input$range[2]) %>%
      group_by_at(dims()) %>%
      summarize(wOBABIP = weighted.mean(wOBABIP, density)) %>%
      list(., n)
  })
  
  wobateam <- reactive({
    req(input$year,
        input$team,
        input$fieldposition)
    calcTeamWoba(input$team,
                 input$awayonly,
                 input$year,
                 input$fieldposition,
                 input$difference)
  })
  
  observe({
    req(input$dim)
    updateSliderInput(
      session,
      "range",
      min = dimRange()[1],
      max = dimRange()[2],
      value = dimRange
    )
  })
  
  output$wobaplotfull <- renderPlot({
    plotWobaMargin(
      wobamarginalized()
    )
  })
  
  output$wobaplotteam <- renderPlot({
    plotWobaTeam(wobateam(),
                 input$difference)
  })
}

shinyApp(ui = ui, server = server)
