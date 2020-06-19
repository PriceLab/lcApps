ID = 18jun2020
IMAGE_NAME = lc_app_$(ID)
TESTPORT=60055

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .

bash:
	docker run -it \
           -p 9044:9044 \
           --entrypoint /bin/bash \
           $(IMAGE_NAME)


testDockerReadiness:
	R -e 'shiny::runApp("lcApp.R")'

