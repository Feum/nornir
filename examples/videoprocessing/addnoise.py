import numpy as np
import os
import cv2
import sys

def noisy(noise_typ, image):
	if noise_typ == "gauss":
		row,col,ch = image.shape
		mean = 0
		var = 0.1
		sigma = var**0.9 # The higher the value the lesser the noise
		gauss = np.random.normal(mean,sigma,(row,col,ch))
		gauss = (gauss.reshape(row, col, ch) * 255.0).astype('u1')
		noisy = image + gauss
		return noisy
	elif noise_typ == "s&p":
		row,col,ch = image.shape
		s_vs_p = 0.5
		amount = 0.004
		out = image
		# Salt mode
		num_salt = np.ceil(amount * image.size * s_vs_p)
		coords = [np.random.randint(0, i - 1, int(num_salt)) for i in image.shape]
		out[coords] = 1

		# Pepper mode
		num_pepper = np.ceil(amount* image.size * (1. - s_vs_p))
		coords = [np.random.randint(0, i - 1, int(num_pepper)) for i in image.shape]
		out[coords] = 0
		return out
	elif noise_typ == "poisson":
		vals = len(np.unique(image))
		vals = 2 ** np.ceil(np.log2(vals))
		noisy = np.random.poisson(image * vals) / float(vals)
		return noisy
	elif noise_typ =="speckle":
		row,col,ch = image.shape
		gauss = np.random.randn(row,col,ch)
		gauss = gauss.reshape(row,col,ch)        
		noisy = image + image * gauss
		return noisy


cap = cv2.VideoCapture(sys.argv[1])
ret, frame = cap.read()
fps = 0
(major_ver, minor_ver, subminor_ver) = (cv2.__version__).split('.')

if int(major_ver)  < 3 :
	fps = cap.get(cv2.cv.CV_CAP_PROP_FPS)
else :
	fps = cap.get(cv2.CAP_PROP_FPS)
     

# Define the codec and create VideoWriter object
fourcc = cv2.cv.CV_FOURCC(*'MPEG')
#out = cv2.VideoWriter('output.avi',fourcc, fps, (640,480))
#fourcc = cap.get(cv2.cv.CV_CAP_PROP_FOURCC) #cv2.cv.CV_FOURCC(*'MJPG')
out = cv2.VideoWriter('output.avi', fourcc, fps, (frame.shape[1], frame.shape[0]))

if not out.isOpened():
	print "Impossible to open the output file."

while(cap.isOpened()):
	ret, frame = cap.read()
	if ret==True:
		frame = noisy("gauss", frame)
		out.write(frame)	
		#cv2.imshow('frame',frame)
		#if cv2.waitKey(1) & 0xFF == ord('q'):
		#	break
	else:
		break

# Release everything if job is finished
cap.release()
out.release()
cv2.destroyAllWindows()
