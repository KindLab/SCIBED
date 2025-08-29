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

# Docker hub

We did this:
```
docker tag scibed-scibed-rstudio robinhweide/scibed-scibed-rstudio:V1
docker login
docker push robinhweide/scibed-scibed-rstudio:V1
```

So you can do:
```
docker pull robinhweide/scibed-scibed-rstudio:V1
```
