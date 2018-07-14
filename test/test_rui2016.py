import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))
from d2d.benchmark import run_rui2016


def main():
    ans = run_rui2016(n_cc=5, n_pairs=5, n_rb=50, cc_qos=1)
    ans.wait()
    print(ans.champion_x)


if __name__ == '__main__':
    main()