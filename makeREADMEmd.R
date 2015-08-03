library(latexreadme)

rmd = "README.Rmd"
md = "README.md"

parse_latex(rmd,md,git_username="IRSN",git_reponame="IECI",insert_string = paste0("\n<img src=\"%s/figure/%s\" height=\"30\">\n"))

for (f in list.files(pattern = "*.png")) {
  #print(f)
  file.copy(f,paste0("figure/",f),overwrite = T)
  file.remove(f)
}

library(knitr)

html = pandoc(md, format = "html")