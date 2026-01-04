library(tidyverse)
library(data.table)
library(corrplot)
library(r.jive)

# Read GGIR files
sleep     <- read_csv("your path", show_col_types = FALSE)
pa        <- read_csv("your path", show_col_types = FALSE)
circadian <- read_csv("your path", show_col_types = FALSE)


# keep eid + all numeric features; drop all-NA and constant columns
prep_domain <- function(df, domain_name, eid_col = "eid") {
  if (!eid_col %in% names(df)) {
    stop("Missing column '", eid_col, "' in ", domain_name, " file.")
  }
  
  # Keep eid + numeric columns only
  out <- df %>%
    dplyr::select(all_of(eid_col), where(is.numeric))
  
  # Drop columns that are all NA (except eid)
  keep1 <- vapply(out, function(x) !all(is.na(x)), logical(1))
  out <- out[, keep1, drop = FALSE]
  
  # Drop constant columns (sd==0 or only one unique finite value), except eid
  num_cols <- setdiff(names(out), eid_col)
  if (length(num_cols) > 0) {
    is_const <- vapply(out[, num_cols, drop = FALSE], function(x) {
      x <- x[is.finite(x)]
      length(unique(x)) <= 1
    }, logical(1))
    
    if (any(is_const)) {
      out <- out %>% dplyr::select(-all_of(names(is_const)[is_const]))
    }
  }
  
  # Basic logging
  cat(domain_name, ": n=", nrow(out), ", p(numeric features)=", ncol(out) - 1, "\n")
  out
}

sleepDomain     <- prep_domain(sleep, "Sleep")
paDomain        <- prep_domain(pa, "Physical Activity")
circadianDomain <- prep_domain(circadian, "Circadian Rhythm")

# Correlation plots (optional)
plot_corr_safe <- function(domain_df, title = "", max_vars = 40) {
  X <- domain_df %>% dplyr::select(-eid)
  if (ncol(X) < 2) {
    message("Skip corrplot (too few numeric vars): ", title)
    return(invisible(NULL))
  }
  if (ncol(X) > max_vars) {
    message(title, " has ", ncol(X), " vars; plotting first ", max_vars, " only.")
    X <- X[, 1:max_vars, drop = FALSE]
  }
  corM <- cor(X, use = "pairwise.complete.obs")
  corrplot(
    corM,
    method = "color",
    type = "lower",
    tl.col = "black",
    tl.srt = 45,
    addCoef.col = "black",
    number.cex = 0.6,
    col = colorRampPalette(c("blue", "white", "red"))(200)
  )
}

plot_corr_safe(sleepDomain, "Sleep domain", max_vars = 40)
plot_corr_safe(paDomain, "PA domain", max_vars = 40)
plot_corr_safe(circadianDomain, "Circadian domain", max_vars = 40)

# Prepare data for JIVE
domainName <- c("Sleep", "Physical Activity", "Circadian Rhythm")

# Merge all domains by eid
jiveData <- sleepDomain %>%
  inner_join(paDomain, by = "eid") %>%
  inner_join(circadianDomain, by = "eid")

# Remove rows with ANY missing values
jiveData_complete <- jiveData %>% drop_na()

# Extract eid for later use
jive.filenames <- jiveData_complete$eid
cat("Complete-case N =", nrow(jiveData_complete), "\n")

# Build JIVE input list: matrices (features x participants)
sleep_vars <- setdiff(names(sleepDomain), "eid")
pa_vars    <- setdiff(names(paDomain), "eid")
cr_vars    <- setdiff(names(circadianDomain), "eid")

jive.data <- list(
  # Sleep domain
  jiveData_complete %>%
    dplyr::select(all_of(sleep_vars)) %>%
    as.matrix() %>%
    t(),
  
  # PA domain
  jiveData_complete %>%
    dplyr::select(all_of(pa_vars)) %>%
    as.matrix() %>%
    t(),
  
  # Circadian domain
  jiveData_complete %>%
    dplyr::select(all_of(cr_vars)) %>%
    as.matrix() %>%
    t()
)

names(jive.data) <- domainName

# Sanity check
cat("JIVE input dims:\n")
for (k in seq_along(jive.data)) {
  cat(" - ", names(jive.data)[k], ": ",
      nrow(jive.data[[k]]), "features x",
      ncol(jive.data[[k]]), "participants\n", sep = "")
}

# Store feature names for each domain
jive.data.rownames <- list(
  colnames(sleepDomain)[-1],
  colnames(paDomain)[-1],
  colnames(circadianDomain)[-1]
)

# Name the list elements
names(jive.data) <- domainName

# Check dimensions
for (k in 1:3) {
  cat(domainName[k], ":", dim(jive.data[[k]]), "\n")
}

# Global transformation function from Kang et al.
global.transform <- function(L0) {
  z <- matrix(0, nrow(L0), ncol(L0))
  
  for (i in 1:ncol(L0)) {
    F <- ecdf(L0[, i])
    u <- F(L0[, i])
    z[, i] <- qnorm(u)
  }
  
  z[is.infinite(z)] = max(z[!is.infinite(z)], na.rm = TRUE) + 1
  colnames(z) <- colnames(L0)
  return(z)
}

# Apply global normalization to each domain
jive.data.norm <- jive.data
for (k in 1:3) {
  jive.data.norm[[k]] <- t(global.transform(t(jive.data.norm[[k]])))
  cat(domainName[k], ":", dim(jive.data.norm[[k]]), "\n")
}

# Run JIVE with automatic rank selection
set.seed(20250601)
res <- jive(jive.data.norm,
            scale = TRUE,
            center = TRUE
)

# View results
summary(res)

# Variance explained plot
showVarExplained(res, col = c("#999999", "#E69F00", "#56B4E9"))

# Heatmap visualization
plot(res, type = "heat")

# Compute JIVE scores
J_Est <- jive.predict(jive.data.norm, res)

joint <- J_Est$joint.scores
rownames(joint) <- paste0("JScore_Joint_", seq_len(nrow(joint)))

indiv_list <- vector("list", 3)
for (k in 1:3) {
  temp <- J_Est$indiv.scores[[k]]
  dn <- gsub("\\s+", "", domainName[k])
  rownames(temp) <- paste0("JScore_", dn, "_", seq_len(nrow(temp)))
  indiv_list[[k]] <- temp
}

all_scores <- rbind(joint, do.call(rbind, indiv_list))

Jscore_df <- data.frame(
  eid = jive.filenames,
  t(all_scores),
  check.names = FALSE
)

# Save scores
write.csv(Jscore_df, file = "jive_predScore.csv", row.names = FALSE)

# PCA plots
showPCA(res, n_joint = 2, Colors = 'purple')
showPCA(
  res,
  n_joint = 1,
  n_indiv = c(1, 1, 1),
  Colors = 'purple'
)

# Save JIVE decomposition
jive.joint <- NULL
jive.ind <- NULL

for (k in 1:3) {
  temp <- cbind(
    jive.data.rownames[[k]],
    domainName[k],
    "J",
    res$joint[[k]]
  )
  jive.joint <- rbind(jive.joint, temp)
  
  temp2 <- cbind(
    jive.data.rownames[[k]],
    domainName[k],
    "A",
    res$individual[[k]]
  )
  jive.ind <- rbind(jive.ind, temp2)
}

jive.out <- rbind(jive.joint, jive.ind)
colnames(jive.out) <- c("feature", "domain", "jive", jive.filenames)

write.csv(jive.out,
          file = "jive_Decomposition.csv",
          row.names = FALSE)
