from utils import *
from scipy.interpolate import griddata
from scipy.signal import medfilt2d, order_filter

def mean_filter(data, dx_dt, window_width):
	
	""" Remove the DC bias by using a mean filter. window_width (ns)"""
	
	dx, dt = dx_dt
	samples, traces = data.shape
	t = np.linspace(0, 1, samples) * (dt * samples)
	
	width = int(window_width / dt)
	half_width = int(width / 2)
	temp = np.ones((samples, traces))
	temp *= data
	
	for trace in range(traces):
		for index in range(samples):
			if index < half_width:
				data[index,trace] += -np.mean(abs(temp[:index + half_width, trace]))
				
			elif index > (samples - half_width):
				data[index,trace] += -np.mean(abs(temp[index-half_width:]))
			else:
				data[index,trace] += -np.mean(abs(temp[index-half_width:index+half_width, trace]))
				
	return data
		
def dc_substraction (data, dx_dt, time_window):
	
	""" Remove the DC bias by using a simple substraction of the DC offset
	time_window is a tuple (start, stop) in ns."""
	dt = dx_dt[1]
	start, stop = time_window
	start = int(start / dt)
	stop = int(stop / dt)
	
	return data - np.mean(abs(data[start:stop]))	
	
def cut_off_frequency(data, dx_dt, fc):
	
	""" Apply a pass-band filter by using Fast Fourier Transform. fc (MHz) must be bellow the bandwidth of the recorded data  ~ 10 MHz """
	
	dx, dt = dx_dt
	samples, traces = data.shape
	f = np.linspace(0, 1, samples) * (samples / dt)
	
	index = int(fc * dt)
	fftdata = np.fft.fft2(data)
	for trace in range(traces):
		fit = np.diff(fftdata.real[index:index+2, trace]) * [range(index)]
		fftdata.real[:index, trace] = fit
	
	data = np.fft.ifft2(fftdata)
	return data.real
	
def time_zero (data, dx_dt, t0 = 0.0):
	
	"""Replaces the start time of your radargrams by t0 (ns), retrun a new 2D dataset reshaped"""
	
	dx, dt = dx_dt
	samples, traces = data.shape
	t = np.linspace(0, 1, samples) * (samples * dt)
	index = int(t0 / dt)

	return data[index:]
		
def user_gain (data, dx_dt, sort, param, time_window, plot=False):
	
	"""Add a user defined gain chosen beetween {'constant', 'linear', 'enxponential'} on data.
	param is a tuple (a, b) where --> constant : fgain = a, linear : fgain = a*t, exponential : fgain = a*exp(b*t)
	time_window is a tuple (start, stop) in ns.
	"""
	dx, dt = dx_dt
	samples, traces = data.shape
	a, b = param
	t = np.linspace(0, 1, samples) * (samples * dt)
	t0, stop = time_window
	
	start = int(t0 / dt)
	stop = int(stop / dt)
	width = start-stop
	fgain = np.ones(samples)
	
	if sort == 'constant':
		fgain[start:stop] = [a]*width
		
	elif sort == 'linear':
		fgain[start:stop] = [a*(t-t0) + 1 for t in t[start:stop]]

	elif sort == 'exponential':
		fgain[start:stop] = [a*(np.exp(b*(t-t0))-1) for t in t[start:stop]]
		
	for trace in range(traces):
		data[:, trace] *= fgain.astype(dtype=data.dtype)
	
	if plot is True:
		plt.plot(fgain)
		plt.show()
	
	return data

def velocity_analysis (data, dx_dt, param, width):
	"""
	Plot the radargram along with the hyperbolic function initialized by a tuple = (x0 (m), t0 (ns), c (m/ns), r (m)).
	"""

	samples, traces = data.shape
	dx, dt = dx_dt
	x0, t0, v, r = param

	mid = int(x0 / dx)
	start = int(mid - width / (2 * dx))
	stop = int(mid + width / (2 * dx))
	
	z0 = (t0 * v + 2 * r) / 2
	x = np.linspace(0, 1, traces) * (traces * dx)
	t = np.linspace(0, 1, samples) * (samples * dt)
	hyperbol =  (2 / v) * (np.sqrt((x0-x[start:stop])**2 + z0**2) - r) 
	
	fig, ax = plt.subplots(num='velocity_analysis', figsize=(20, 10), facecolor='w', edgecolor='w')
	ax.imshow(data, extent=[np.amin(x), np.amax(x), np.amax(t), np.amin(t)], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(abs(data)), vmax=np.amax(abs(data)))
	ax.plot(x[start:stop], hyperbol)
	plt.show()
	
def stolt_migration(data, dx_dt, c):
	
	# Imput parameters
	dx, dt = dx_dt
	fs = 1/dt
	eps = 2.2e-16
	nt0, nx0 = data.shape
	t = np.linspace(0,nt0*dt,nt0) 
	x = np.linspace(0,nx0*dx,nx0)
	
	# Zero-padding 2D
	nt = 2 * nextpower(nt0)
	nx = 2 * nx0
	
	# One Emiter-Receiver scenario
	ERMv = c / 2
	
	# FFT & shift 
	fftdata = np.fft.fftshift(np.fft.fft2(data, s=(nt,nx)))
	
	# Build (kx, f) 
	f = np.linspace(-nt/2, nt/2-1, nt) * fs / nt 
	kx = np.linspace(-nx/2,nx/2-1, nx) / dx / nx
	kx, f = np.meshgrid(kx, f)
	
	# Remove evanescent parts
	evanescent = (abs(f)  / abs(kx+eps) > c).astype(int)
	fftdata *= evanescent
	
	# Stolt remapping function f(kz)
	fkz = ERMv*np.sign(f)*np.sqrt(kx**2 + f**2/ERMv**2)
	
	# Linear interpolation on grid
	fftdata = griddata((kx.ravel(), f.ravel()), fftdata.ravel(), (kx, fkz), method='nearest')
	
	# Jacombien
	fftdata *= f / (np.sqrt(kx**2 + f**2/ERMv**2)+eps)
	
	# IFFT & Migrated data 
	mig = np.fft.ifft2(np.fft.ifftshift(fftdata))
	mig = mig[:nt0,:nx0]
	
	dz = dt * c / 2
	dx_dz = (dx, dz)
	
	return abs(mig), dx_dz
	
	
	
	
