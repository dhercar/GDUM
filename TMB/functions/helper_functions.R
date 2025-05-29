print_title <- function(text, width = getOption("width"), symb = '-') {
  # Full-width line of dashes
  separator_line <- strrep(symb, width)
  
  # Centered text
  padding <- floor((width - nchar(text)) / 2)
  if (padding < 0) padding <- 0
  centered_text <- paste0(strrep(" ", padding), text)
  
  # Print all
  cat(separator_line, "\n")
  cat(centered_text, "\n")
  cat(separator_line, "\n\n")
}

print_title2 <- function(text, width = getOption("width"), symb = '-') {
  total_padding <- width - nchar(text)
  if (total_padding < 0) total_padding <- 0
  
  left_padding <- floor(total_padding / 2)
  right_padding <- ceiling(total_padding / 2)
  
  decorated_text <- paste0(strrep(symb, left_padding), text, strrep(symb, right_padding))
  cat(decorated_text, "\n\n")
}


