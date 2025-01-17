test_that("Test design construction", {

  events_df <- data.frame(
    onsets = c(8.8, 24.2, 55, 70.4, 85.8, 101.2, 132, 147.4, 162.8, 193.6),
    duration = rep(15.4, 10),
    trial_type = c("RHand", "RFoot", "LFoot", "Tongue", "LHand", "LFoot", "LHand", "RHand", "Tongue", "RFoot")
  )

  ttypes <- sort(unique(events_df$trial_type))
  events <- setNames(lapply(
    ttypes,
    function(ttype) { events_df[events_df$trial_type==ttype,c("onsets", "duration")] }
  ), ttypes)

  # Make version of `events` with no stimuli for some tasks, for testing.
  eventsB <- events
  eventsB[c("LHand", "RHand")] <- NA

  nTime <- 350
  TR <- .72

  des <- make_design(events, nTime, TR)

  testthat::expect_no_error(
    print(des)
  )

  testthat::expect_no_error(
    plot(des)
  )

  testthat::expect_warning(
    desB <- make_design(eventsB, nTime, TR)
  )
})
