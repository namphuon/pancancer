 # Use an official Python runtime as a parent image
FROM ubuntu

# Set the working directory to /app
WORKDIR /home/
RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
#Download libraries for AA
RUN apt-get update && apt-get install -y
RUN apt-get install python-dev gfortran python-numpy python-scipy python-matplotlib python-pip zlib1g-dev -y
RUN pip install pysam gcloud

#Set up gcloud
RUN apt-get install lsb-release curl -y
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)"
#RUN echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
RUN echo "deb http://packages.cloud.google.com/apt cloud-sdk-stretch main" >  /etc/apt/sources.list.d/google-cloud-sdk.list
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
RUN apt-get update && apt-get install google-cloud-sdk -y
RUN apt-get install gcc python-dev python-setuptools libbz2-dev liblzma-dev -y
RUN easy_install -U pip
RUN pip install -U crcmod
RUN pip install --upgrade google-api-python-client
RUN pip install requests
ENV PANCANCER=/home

ADD *.py /home/
ADD read_count.sh /home/
ADD *.bz2 /home/
RUN cd /home/htslib-1.7 && ./configure && make && make install
RUN cd /home/bcftools-1.7 && ./configure && make && make install
RUN apt-get install libncurses5-dev libbz2-dev liblzma-dev -y
RUN cd /home/samtools-1.7 && ./configure && make && make install