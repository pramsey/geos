# GEOS Web Site

Built using [Jekyll](https://jekyllrb.com/).

To run locally while editing, the easiest path is to use the docker image, and the following command:

```
export JEKYLL_VERSION=3.8
docker run --rm \
  --volume="$PWD:/srv/jekyll" -p 4000:4000 \
  --volume="$PWD/vendor/bundle:/usr/local/bundle" \
  -it jekyll/builder:$JEKYLL_VERSION \
  jekyll serve \
    --incremental \
    --drafts
```

