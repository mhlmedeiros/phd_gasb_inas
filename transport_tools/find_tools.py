
import numpy as np

def findExtr(xArray, yArray, indIni, delSign = 1):

    """
    Find the extreme value location of a discrete function given by a map
    between xArray and yArray.

    xArray  = array with the domain points of the function
    yArray  = array with the image points
    indIni  = first value index: where to start the search
    delSign = variation signal : 1 -> Maximum, 0 -> Minimum

    """

    for i in range(indIni, len(xArray) - 1):
        delta = np.sign(yArray[i+1] - yArray[i])

        if delta != delSign:
            # Treat the special case when the difference is zero: simmetric points
            if (delta == 0) and (i + 2 < len(xArray)):
                delta = sign(yArray[i+2] - yArray[i+1])
            return i, delta

    return len(xArray), delta

def findAllExtr(xArray, yArray, delSignIni = 1):

    extrIndexes = []
    indIni = 0

    while indIni < len(xArray):
        # Find point where the sign of the differences changes
        indIni, delSignIni = findExtr(xArray, yArray, indIni, delSignIni)
        # Append into the list
        extrIndexes.append(indIni)

    return extrIndexes

def main():
    xname = "data_448_450.1_meV_Fermi_2001_L_707.1067811865476"
    yname = "data_448_450.1_meV_Fermi_Transport_Total_2001_L_707.1067811865476"
    xArray = np.load("../data/transport/" + xname + ".npy")
    yArray = np.load("../data/transport/" + yname + ".npy")

    reverso = False

    if reverso == True:
        xArray = xArray[::-1]
        yArray = yArray[::-1]

    x_lim_min = 448
    x_lim_max = 450
    

    x_new = xArray[(xArray >= x_lim_min) & (xArray <= x_lim_max)]
    y_new = yArray[(xArray >= x_lim_min) & (xArray <= x_lim_max)]

    indexes = findAllExtr(x_new, y_new)

    # print(indexes)
    xExtrema = np.asarray([xArray[j] for j in indexes[:-1]])
    yExtrema = np.asarray([yArray[j] for j in indexes[:-1]])

    np.save("../data/transport/" + 
                xname + "_extrema_Fermi_" +
                str(x_lim_min) + "_" + str(x_lim_max) + ".npy", xExtrema)
    np.savetxt("../data/transport/" + 
                xname + "_extrema_Fermi_" + 
                str(x_lim_min) + "_" + str(x_lim_max) + ".txt", xExtrema)

    np.save("../data/transport/"+ yname + "_extrema_Fermi_" +
                str(x_lim_min) + "_" + str(x_lim_max) + ".npy", yExtrema)
    np.savetxt("../data/transport/"+ yname + "_extrema_Fermi_" +
                str(x_lim_min) + "_" + str(x_lim_max) + ".txt", yExtrema)

if __name__ == "__main__":
    main()
