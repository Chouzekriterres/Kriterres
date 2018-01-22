import theta2eps
import numpy as np
import matplotlib.pyplot as plt


def main():
    print("Hello World!")
    theta = np.arange(0.06, 0.4, 0.05)
    print(theta)
    porosity = 0.45
    eps = theta2eps.CRIM(theta, porosity)
    plt.plot(theta, eps, 'g')
    plt.show()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\epsilon$')


if __name__ == "__main__":
    main()
