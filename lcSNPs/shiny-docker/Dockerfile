FROM bioconductor/bioconductor_docker
#FROM rocker/shiny
WORKDIR /home/shiny
COPY dataTableModule.R .
COPY iframeSearchModule.R .
COPY igvModule.R .
COPY installPrerequisites.R .
COPY lcApp.R .
COPY messageBoxModule.R .
COPY tbl.gh.apoe.RData .
COPY threeSnps.tsv .
RUN R -f ./installPrerequisites.R .
COPY pkgs/VariantAnnotation_1.34.0.tar.gz .
RUN R CMD INSTALL VariantAnnotation_1.34.0.tar.gz
COPY pkgs/igvShiny_0.99.90.tar.gz .
RUN R CMD INSTALL igvShiny_0.99.90.tar.gz
CMD ["R", "-e", "shiny::runApp('/home/shiny/lcApp.R', port=3838, host='0.0.0.0')"]





