FROM lc_app_base_19jun2020
WORKDIR /home/shiny
COPY lcApp.R .
COPY tbl.gh.apoe.RData .
COPY tbl.snps.598.RData .
COPY tbl.snps.630.RData .
COPY tbl.summary.1010x6-gene-rowNames.RData .
COPY pkgs/igvShiny_0.99.90.tar.gz .
RUN R CMD INSTALL igvShiny_0.99.90.tar.gz
COPY pkgs/shinyModules_1.0.4.tar.gz .
RUN R CMD INSTALL shinyModules_1.0.4.tar.gz
# CMD ["R", "-e", "shiny::runApp('/home/shiny/lcApp.R', port=3838, host='0.0.0.0')"]





