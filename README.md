# Hotspots
predict hotspots for MHC-I ligands in proteome

## pipeline description ##
* apply sliding window over the proteome and extract the mean count over a window as the window's count
* split the data into train and test data set
* calculate one-hot and normalised AA index encodings for each window
* fit a deep convolutional neural net to the training data and run a prediction on the test data
* evaluate the performance (PCC, R^2, profile plots, ...)
