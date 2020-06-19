FROM lc_app_base_19jun2020
WORKDIR /home/shiny
COPY dataTableModule.R .
COPY iframeSearchModule.R .
COPY igvModule.R .
COPY lcApp.R .
COPY messageBoxModule.R .
COPY tbl.gh.apoe.RData .
COPY threeSnps.tsv .
COPY pkgs/igvShiny_0.99.90.tar.gz .
RUN R CMD INSTALL igvShiny_0.99.90.tar.gz
# CMD ["R", "-e", "shiny::runApp('/home/shiny/lcApp.R', port=3838, host='0.0.0.0')"]





