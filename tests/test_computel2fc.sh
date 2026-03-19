#### TEST SUITE ####
#
#  Use build_docker.sh to build the docker image.
#  Make sure to change the docker tag to "develop"
#  Please minimize assets before adding to repo, and delete output after testing
#
#####
cd ../
docker build -t prismtheprismlab/sushi:develop --rm=true -f docker/Dockerfile .

DOCKER_TAG="develop"

docker run -it -v /Users/anup/work/:/Users/anup/work/ prismtheprismlab/sushi:$DOCKER_TAG compute_l2fc \
--normalized_counts /Users/anup/work/1_Projects/AWS_sushi/EPS-001_reprocessed2/normalized_counts.csv \
--count_col_name normalized_n --control_type negcon --sig_cols cell_set,treatment,dose,dose_unit,day \
--ctrl_cols cell_set,day --out /Users/anup/work/1_Projects/AWS_sushi/EPS-001_reprocessed2/

