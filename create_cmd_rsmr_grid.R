create_cmd_rsmr_grid <- function(cmd) {
  # Filter by PERIOD
  CMD <- cmd 
  #%>% select(BGC, all_of(cmd.type))
  
  # Compute wetter rSMR values
  for (i in 1:3) {
    CMD[, 2 + i] <- CMD[, 1 + i] / 2
  }
  
  # Rename columns
  colnames(CMD) <- c("BGC", "rSMR4", "rSMR5", "rSMR6", "rSMR7")
  CMD <- CMD[, c("BGC", "rSMR7", "rSMR6", "rSMR5", "rSMR4")]
  # Compute drier rSMR values
  for (i in 1:4) {
    CMD[, length(CMD) + 1] <- CMD[, length(CMD)] + 125
  }
  
  # Rename newly added columns
  colnames(CMD)[6:9] <- c("rSMR3", "rSMR2", "rSMR1", "rSMR0")
  
  # Reorder columns
  CMD <- CMD[, c("BGC", "rSMR7", "rSMR6", "rSMR5", "rSMR4", "rSMR3", "rSMR2", "rSMR1", "rSMR0")]
  
  # Filter by specified BGCs
  #CMD <- CMD %>% filter(BGC %in% AllBGC)
  
  # Add cmd.type and period identifiers
  #CMD$cmd.type <- cmd.type
  #CMD$PERIOD <- period
  
  return(CMD)
}