library(simmer)
set.seed(999)

env <- simmer("mySim")

patient <- trajectory("patients' path") %>%
seize("nurse", 1) %>%
timeout(function() rnorm(1, 15)) %>%
release("nurse", 1) %>%
seize("doctor", 1) %>%
timeout(function() rnorm(1, 20)) %>%
release("doctor", 1) %>%
seize("administration", 1) %>%
timeout(function() rnorm(1, 5)) %>%
release("administration", 1)

env %>%
add_resource("nurse", 1) %>%
add_resource("doctor", 2) %>%
add_resource("administration", 1) %>%
add_generator("patient", patient, function() rnorm(1, 10, 2))

env %>%
  run(80) %>%
  now()