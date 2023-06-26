ID = 19may2020
IMAGE_NAME = lc_genes_$(ID)
TESTPORT=60055

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .

bash:
	docker run -it \
           -p 9012:3838 \
           --entrypoint /bin/bash \
           $(IMAGE_NAME)


shinyRun:
	docker run -it \
           -v ~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos/shiny/data:/home/shiny/data \
           -p 9010:3838 \
           $(IMAGE_NAME)
	docker run -ti --rm -p $(TESTPORT):3838 rocker/shiny bash

dockerRun:
	docker run -p $(TESTPORT):3838 $(IMAGE_NAME)



testDockerReadiness:
	R -e 'shiny::runApp("lcGenes.R")'

