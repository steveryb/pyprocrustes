from scipy import stats
import numpy as np

import csv
data_sets = []
filenames = ["q_3d.csv","q_s.csv"]
for filename in filenames:
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        data = []
        for row in reader:
            data.append(row)
        data = np.transpose([[int(x) if x else 0 for x in row[1:]] for row in data[1:]])
        data_sets.append([[x for x in row if x!=0] for row in data])
print("mean 3d,mean sky,diff mean,t-test Welch t-test, Welch two-tailed prob")
for row in range(len(data_sets[0])):
    d_row = data_sets[0][row]
    sky_row = data_sets[1][row]
    d_mean = np.mean(d_row)
    sky_mean = np.mean(sky_row)
    diff_mean = d_mean-sky_mean
    t_test = stats.ttest_ind(d_row,sky_row)
    w_t_test = stats.ttest_ind(d_row,sky_row,equal_var=False)
    print row+1,d_mean,sky_mean,diff_mean,w_t_test[0],w_t_test[1]
