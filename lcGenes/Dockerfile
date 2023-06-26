FROM rocker/shiny
WORKDIR /home/shiny
COPY tbl.summary.1010x6.RData .
COPY tbl.orthologsBySpecies.9999x4.RData .
COPY geneTableHelp.html .
COPY introduction.html .
COPY installPrerequisites.R .
COPY lcGenes.R .
COPY lcGenes.css .
RUN R -f ./installPrerequisites.R .
CMD ["R", "-e", "shiny::runApp('/home/shiny/lcGenes.R', port=3838, host='0.0.0.0')"]





