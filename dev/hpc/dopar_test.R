# Quick test that foreach + doParallel actually parallelizes on HPC.
# Each iteration sleeps 1s. With 4 workers running 8 jobs sequentially
# you'd see ~8s; with parallelism, ~2s.

library(foreach)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

shared_var <- "hello from the main process"

t0 <- Sys.time()
results <- foreach(i = 1:8, .combine = rbind) %dopar% {
  Sys.sleep(1)
  data.frame(
    iter   = i,
    pid    = Sys.getpid(),
    sees   = shared_var,
    stamp  = format(Sys.time(), "%H:%M:%OS3")
  )
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

stopCluster(cl)

cat("\n--- Results ---\n")
print(results)
cat("\nUnique worker PIDs: ", length(unique(results$pid)),
    " (should be > 1 if parallel)\n", sep = "")
cat("Elapsed: ", round(elapsed, 1), "s",
    " (sequential = 8s; parallel w/ 4 workers = ~2s)\n", sep = "")

if (elapsed > 5) {
  cat("\n[!] Looks SEQUENTIAL — foreach is not actually parallelizing.\n")
} else {
  cat("\n[OK] Parallel — foreach + doParallel works on this node.\n")
}
