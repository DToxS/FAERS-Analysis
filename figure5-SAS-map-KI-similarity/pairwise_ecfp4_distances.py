###generates pairwise tanimoto coefficients given an input smiles file (formatted: SMILES,lig_id)
from rdkit import DataStructs
from rdkit import Chem 
import rdkit 
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors
import sys
smi = Chem.SmilesMolSupplier(sys.argv[1],delimiter=',',titleLine=True)
fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, useBondTypes=False , nBits=1024) for x in smi] ### ECFP4
fps2 = [AllChem.GetMorganFingerprintAsBitVect(x,1, useBondTypes=False , nBits=1024) for x in smi] ### ECFP2 
maccs = [MACCSkeys.GenMACCSKeys(x) for x in smi]
dl = [rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(x) for x in smi]
print("D1,D2,ECFP4,ECFP2,MACCS,DL,AVG,Weighted")
seen = []
for i in range(len(fps)): 
    d1 = smi[i].GetProp('_Name')
    for ii in range(len(fps)):
        d2 = smi[ii].GetProp('_Name')
        dist  = DataStructs.FingerprintSimilarity(fps[i], fps[ii])
        dist2  = DataStructs.FingerprintSimilarity(fps2[i], fps2[ii])
        distMACCS  = DataStructs.FingerprintSimilarity(maccs[i], maccs[ii])
        distDL  = DataStructs.FingerprintSimilarity(dl[i], dl[ii])
        weightedavg = dist*.3 + dist2*.3 + distDL*.3 + distMACCS*.1 
        avg = ( dist + dist2 + distDL + distMACCS ) / 4
        print(str(d1) + "," + str(d2) + "," +  str(dist) + "," +  str(dist2) + "," + str(distMACCS) + ","  + str(distDL) + ","  + str(avg) + ","  + str(weightedavg))
    seen.append(d1)
