ID = 23jun2020
BASE_IMAGE_NAME = lc_app_base_$(ID)
APP_IMAGE_NAME = lc_app_$(ID)

build:
	docker build -t $(APP_IMAGE_NAME) -f Dockerfile .

build-base:
	docker build -t $(BASE_IMAGE_NAME) -f Dockerfile-base .

bash:
	docker run -it \
           -p 9051:3838 \
           --entrypoint /bin/bash \
           $(APP_IMAGE_NAME)

