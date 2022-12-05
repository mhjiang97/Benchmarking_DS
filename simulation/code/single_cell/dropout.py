import os, math
import numpy as np
import pandas as pd

def logistic(x):
    return np.exp(x) / (1 + np.exp(x))

def logit(x, minval = 0.001):
    if isinstance(x, (list, tuple, np.ndarray)):
        x[1 - x < minval] = 1 - minval
        x[x < minval] = minval
    else:
        x = max(minval, x)
        x = min(1 - minval, x)
    val = np.log(x / (1 - x))
    return val

def adjust_drop_prob(drop_prob, rate_new = 0.3):
    gaps_all = np.arange(-10, 10, 0.05)
    rate_all = np.zeros(len(gaps_all))
    drop_logit = logit(drop_prob)
    for i in range(len(gaps_all)):
        drop_prob_tmp = logistic(drop_logit + gaps_all[i])
        rate_all[i] = np.mean(drop_prob_tmp)
    idx = np.argmin(np.abs(rate_all - rate_new))
    drop_prob_new = logistic(drop_logit + gaps_all[idx])
    return drop_prob_new


def fpkmToTPM(fpkm):
    total = sum(fpkm)
    tpm = [math.exp(math.log(x) - math.log(total) + math.log(1e6)) if x != 0 else 0 for x in fpkm]
    return tpm

def dropOut(data_frame, probability = 0.001, dropout_rate = None):
    df = data_frame.copy()
    fpkm_all = np.array(df.FPKM)
    idx_drop = fpkm_all > 0
    fpkm_fil = fpkm_all[idx_drop]
    dropout_prob = np.ones(len(fpkm_fil)) * probability
    if dropout_rate is not None:
        dropout_prob = adjust_drop_prob(dropout_prob, dropout_rate)
    fpkm_do = np.zeros(len(fpkm_fil))
    for i in range(len(fpkm_fil)):
        keep = np.random.binomial(1, 1 - dropout_prob[i])
        fpkm_do[i] = round(keep * fpkm_fil[i], 2)
    df.loc[df.FPKM > 0, "FPKM"] = fpkm_do
    tpm = fpkmToTPM(df.FPKM)
    df["TPM"] = [round(x, 2) for x in tpm]
    df["FPKM"] = [str(x) for x in df.FPKM]
    df["TPM"] = [str(x) for x in df.TPM]
    df["IsoPct"] = [str(x) for x in df.IsoPct]
    df["expected_count"] = [str(x) for x in df.expected_count]
    df["effective_length"] = [str(x) for x in df.effective_length]
    df.loc[df.FPKM == "0.0", "FPKM"] = "0.00"
    df.loc[df.TPM == "0.0", "TPM"] = "0.00"
    df.loc[df.IsoPct == "0.0", "IsoPct"] = "0.00"
    df.loc[df.expected_count == "0.0", "expected_count"] = "0.00"
    df.loc[df.effective_length == "0.0", "effective_length"] = "0.00"
    return df


def main():
    dir_dat = "../../profiles/single_cell/100/"
    for i in range(1, 41):
        # if i == 2:
        #     break
        s = "sample{}".format(str(i)) # s = "sample1"
        print("Processing {}: {}/{}.txt.\n".format(s, dir_dat, s))
        mydf = pd.read_table(os.path.expanduser("{}/{}.txt".format(dir_dat, s)), sep = "\t", header = 0)
        mydf_do = dropOut(data_frame = mydf, dropout_rate = 0.08) # , probability = 0.05
        mydf_do.to_csv(os.path.expanduser("{}/dropout/{}.txt".format(dir_dat, s)), sep = "\t", header = True, index = False, mode = "w")
        print("Done.\n")

if __name__ == "__main__":
    main()




