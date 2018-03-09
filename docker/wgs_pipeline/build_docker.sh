docker build . -t wgs_pipeline
docker tag read_count us.gcr.io/aa-test-175718/wgs_pipeline
gcloud docker -- push us.gcr.io/aa-test-175718/wgs_pipeline