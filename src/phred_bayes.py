fr


if __name__ == '__main__':



single_site_results = [('A',10),('A',5),('G',11),('C',2)]
single_site_results = [('C',2),('A',10),('A',5),('G',11)]
single_site_results = [('G',3),('C',3),('A',3),('A',3),('A',3)]

index_dict = {'A':0,'G':1,'C':2,'T':3}

dist = init_prior()

for read in single_site_results:
    dist = update_dist(dist, get_likelihood(read))
    print(dist)
    print(-np.log(dist)/np.log(10))