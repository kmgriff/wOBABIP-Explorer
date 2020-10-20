##########################################################
# wOBABIP.R
# This script is used for initially pulling and updating raw data, and for precomputing
# large regression objects.  Unless these files need to be updated, you can proceed directly
# to bipexplorer.R to generate plots.
#
# Performs a nonparametric kernel regression on MLB Statcast data with launch angle, launch
# speed, and spray angle as the random variables, with wOBA (weighted on-base average) as the
# outcome variable.  Statcast provides data going back to about 2017.
#
# Bandwidth calculations for the kernel regression and unconditional density estimator, as
# well as the regression and density calculation themselves, all run in O(n^2) time.  The
# bandwidths are calculated on a small random sample of the data, n = 3000, while the
# regression is performed on a grid of specified size corresponding to the density of actual
# data points.  This approach has been roughly approved by the author of the 'np' package as
# statistically sound.  The correct bandwidths computed for this data are roughly 5 for
# vertical launch angle, 50 for launch speed, and 500 for horizontal launch angle.
##########################################################


library(tidyverse)
library(lubridate)
library(baseballr)
library(np)
library(sp)
library(geometry)

ftime <- function(f) {
  f <- as_mapper(f)
  function(...) {
    pre.time <- proc.time()
    x <- f(...)
    print(proc.time() - pre.time)
    x
  }
}

# Collect Statcast data

# Pull initial Statcast data in chunks of 3 days, store in current working directory as "all_mlb.RDS"
if (!"all_mlb.RDS" %in% dir("data")) {
  rawdat <- data.frame()
  for (dates in seq(as_date("2017-01-01"), Sys.Date(), 3)) {
    tdat <- scrape_statcast_savant(dates, dates + 2) %>%
      mutate_at("game_type", as.character)
    if (length(tdat[[1]]) != 0) {
      rawdat <- bind_rows(rawdat, tdat)
    }
  }
  
  # Read data from "all_mlb.RDS", remove last day of data, then update
} else {
  rawdat <- readRDS("data/all_mlb.RDS") %>%
    filter(game_date < Sys.Date() - 1)
  for (dates in seq(max(rawdat$game_date), Sys.Date(), by = "days")) {
    tdat <-
      as_date(dates) %>%
      scrape_statcast_savant(., .) %>%
      mutate_at("game_type", as.character)
    if (length(tdat[[1]]) != 0) {
      rawdat <- bind_rows(rawdat, tdat)
    }
  }
}
rm(tdat, dates)
rawdat <- rawdat %>%
  filter(type == "X") %>%
  filter_at(vars("launch_speed", "launch_angle", "hc_x", "hc_y"),
            ~ !is.na(.x))


rawdat <- rawdat %>%
  mutate(spray_angle = (-atan((hc_x - 130) / (213 - hc_y)) + pi / 2) / pi * 180) %>%
  mutate(spray_angle = 180 - spray_angle) %>%
  filter_at(vars("spray_angle", "woba_value"), ~ !is.na(.x) &
              !is.nan(.x)) %>%
  mutate(fielding_team =
           case_when(
             inning_topbot == "Top" ~ home_team,
             inning_topbot == "Bot" ~ away_team
           )) %>%
  mutate(hitting_team =
           case_when(
             inning_topbot == "Top" ~ away_team,
             inning_topbot == "Bot" ~ home_team
           ))

saveRDS(rawdat, "data/all_mlb.RDS")

# Used for exploration only
calcWobaBandwidths <- ftime(function(n) {
  print(paste0("Calculating wOBA Bandwidths 3D: ", n))
  rawdat %>%
    slice_sample(n = n) %>%
    {
      npregbw(
        xdat = data.frame(
          x = .$launch_angle,
          y = .$launch_speed,
          z = .$spray_angle
        ),
        ydat = as.numeric(.$woba_value),
        bws = c(5, 50, min(n - 1, 500)),
        nmulti = 1,
        bwtype = "generalized_nn",
        bwmethod = "cv.aic"
      )
    }
})

# Used for exploration only
calcWobaBandwidths2d <- ftime(function(n) {
  print(paste0("Calculating wOBA Bandwidths 2D: ", n))
  rawdat %>%
    slice_sample(n = n) %>%
    {
      npregbw(
        xdat = data.frame(x = .$launch_angle, y = .$launch_speed),
        ydat = as.numeric(.$woba_value),
        # bws = c(5, 50),
        # nmulti = 1,
        bwtype = "adaptive_nn",
        bwmethod = "cv.aic"
      )
    }
})

# Used for exploration only
calcDensityBandwidths <- ftime(function(n) {
  print(paste0("Calculating Density Bandwidths 3D: ", n))
  rawdat %>%
    slice_sample(n = n) %>%
    {
      npudensbw(
        dat = data.frame(
          x = .$launch_angle,
          y = .$launch_speed,
          z = .$spray_angle
        ),
        bws = c(5, 50, min(n - 1, 500)),
        nmulti = 1,
        bwtype = "generalized_nn",
        bwmethod = "cv.ml"
      )
    }
})

# Used for exploration only
calcDensityBandwidths2d <- ftime(function(n) {
  print(paste0("Calculating Density Bandwidths 2D: ", n))
  rawdat %>%
    slice_sample(n = n) %>%
    {
      npudensbw(
        dat = data.frame(x = .$launch_angle, y = .$launch_speed),
        bws = c(5, 50),
        nmulti = 1,
        bwtype = "generalized_nn",
        bwmethod = "cv.ml"
      )
    }
})

# bws2d <- calcWobaBandwidths2d(3000)
# bws3d <- calcWobaBandwidths(3000)
# densbws2d <- calcDensityBandwidths2d(3000)
# densbws3d <- calcDensityBandwidths(3000)

dat <- rawdat

# Parse data into a defined grid
grid <- expand.grid(
  launch_angle = seq(-90, 90, 2),
  launch_speed = seq(20, 140, 2),
  spray_angle = seq(45, 135, by = 10)
)
grid2d <-
  expand.grid(launch_angle = seq(-90, 90, 2),
              launch_speed = seq(20, 140, 2))

txdat <-
  data.frame(
    x = dat$launch_angle,
    y = dat$launch_speed,
    z = dat$spray_angle
  )

# Remove grid coordinates that fall outside of convex hull
hull <- convhulln(as.matrix(txdat))
grid_inhull <- grid[inhulln(hull, as.matrix(grid)),]
saveRDS(grid_inhull, "data/gridfull.RDS")

# Calculate mean wOBA at each grid coordinate
calcGridWobaFull <- ftime(function(grid) {
  print("Calculating mean wOBA for grid")
  npreg(
    bws = c(5, 50, 500),
    bwtype = "generalized_nn",
    exdat = grid,
    txdat = txdat,
    tydat = as.numeric(dat$woba_value)
  )
})
gridwobafull <- calcGridWobaFull(grid_inhull)
out <- list(gridwobafull$mean, gridwobafull$ntrain)
names(out) <- c("mean", "n")
saveRDS(out, "data/gridwobafull.RDS")



# Calculate density of observations at each grid coordinate
calcGridDensFull <- ftime(function(grid) {
  print("Calculating density for grid")
  grid %>%
    map( ~ {
      attributes(.x) <- NULL
      .x
    }) %>%
    as.data.frame() %>%
    npudens(
      bws = c(5, 40, 500),
      bwtype = "generalized_nn",
      edat = .,
      tdat = txdat
    )
})
griddens <- calcGridDensFull(grid_inhull)
out <- list(griddens$dens, griddens$ntrain)
names(out) <- c("dens", "n")
saveRDS(out, "data/griddensfull.RDS")



calcWobaDens2d <- function(team, position, away, year) {
  subset <- dat %>%
    filter(lubridate::year(game_date) == year)
  if(position == "hitting") {
    subset <- subset %>% filter(hitting_team == team)
  } else {
    subset <- subset %>% filter(fielding_team == team)
  }
  
  if(away == "away") {
    subset <- subset %>% filter(away_team == team)
  }
  
  txdat2d <- data.frame(x = subset$launch_angle,
                        y = subset$launch_speed)
  
  # Remove grid coordinates that fall outside of convex hull
  hull2d <- convhulln(as.matrix(txdat2d))
  grid_inhull2d <- grid2d[inhulln(hull2d, as.matrix(grid2d)),]
  
  # Calculate mean wOBA at each grid coordinate
  print(str_glue("Calculating mean wOBA 2D for {team} {position} in {away} games in {year}"))
  gridwoba2d <- ftime(function() npreg(
    bws = c(5, 40),
    bwtype = "adaptive_nn",
    exdat = grid_inhull2d,
    txdat = txdat2d,
    tydat = as.numeric(subset$woba_value)
  ))()
  
  # Calculate density of observations at each grid coordinate
  print(str_glue("Calculating density 2D for {team} {position} in {away} games in {year}"))
  griddens <- ftime(function() grid_inhull2d %>%
                      map( ~ {
                        attributes(.x) <- NULL
                        .x
                      }) %>%
                      as.data.frame() %>%
                      npudens(
                        bws = c(5, 40),
                        bwtype = "adaptive_nn",
                        edat = .,
                        tdat = txdat2d
                      ))()
  out <- list(bind_cols(grid_inhull2d,
                        mean = gridwoba2d$mean,
                        dens = griddens$dens),
              gridwoba2d$ntrain)
  names(out) <- c("reg", "n")
  saveRDS(out, str_glue("data/reg2d_{team}_{position}_{away}_{year}.RDS"))
}

walk(dat %>% distinct(year(game_date)) %>% pull(), function(year){
  walk(dat$hitting_team %>% unique() %>% sort(), function(team) {
    walk(c("away", "all"), function(away) {
      walk(c("hitting", "fielding"), function(position) {
        calcWobaDens2d(team, position, away, year)
      })
    })
  })
})
