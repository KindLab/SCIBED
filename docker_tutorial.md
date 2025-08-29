# Usage from scratch

1. Clone your repo (for docker-compose file):
```
git clone https://github.com/KindLab/SCIBED.git
cd SCIBED
```

2. Build and run:
`docker-compose up --build`
Open browser: http://localhost:8787
```
Login:
  user: rstudio
  password: scibed
```

# Usage from image

Since we do `docker save -o scibed_docker_image.tar scibed-scibed-rstudio` with every update, you can load the image directly like so:

```
docker load -i scibed_docker_image.tar
```
