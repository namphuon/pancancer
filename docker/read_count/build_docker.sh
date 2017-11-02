cp ../../python/read_count.py .
cp ../../python/pancancer.py .
docker build . -t read_count
docker tag read_count us.gcr.io/aa-test-175718/read_count
gcloud docker -- push us.gcr.io/aa-test-175718/read_count
