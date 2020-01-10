
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
    xname = "data_62_64_eF_501_clean"
    yname = "data_62_64_transport_501_clean"
    xArray = np.load("../data/transport/" + xname + ".npy")
    yArray = np.load("../data/transport/" + yname + ".npy")


    indexes = findAllExtr(xArray, yArray)

    print(indexes)
    xExtrema = np.asarray([xArray[j] for j in indexes[:-1]])
    yExtrema = np.asarray([yArray[j] for j in indexes[:-1]])

    np.save("../data/transport/"+ xname + "_extrema.npy", xExtrema)
    np.savetxt("../data/transport/"+ xname + "_extrema.txt", xExtrema)

    np.save("../data/transport/"+ yname + "_extrema.npy", yExtrema)
    np.savetxt("../data/transport/"+ yname + "_extrema.txt", yExtrema)

if __name__ == "__main__":
    main()
