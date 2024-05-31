import ecc_solution as ecc

# import A2_ecc_template as ecc
import random
import traceback

print("[*] Teste ECC-Klassen")

try:
    curve = ecc.Curve(
        0xA9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377,
        0x7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9,
        0x26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6,
        0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,
        0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997,
        0xA9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7,
    )
except Exception:
    print(f"[-] Curve konnte nicht erstellt werden")
    traceback.print_exc()
    print("[!] Ende")
    exit(-1)

scalar = 0xFAFAFDEADBEEFAFAFAFAFAFAFAF
try:
    point = ecc.Point(
        0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,
        0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997,
        curve,
    )
except Exception:
    print(f"[-] Point konnte nicht erstellt werden")
    traceback.print_exc()
    print("[!] Ende")
    exit(-1)

if not point.oncurve():
    print("oncurve() liefert fehlerhafte Ergebnisse")
    print("[!] Ende")
    exit(-1)

compare_val = ecc.Point(
    0x353695611afdda8138dff99456ca7c2c228119f82f43cb39ca6fe206e6bde384,
    0x3b87ca7460ea1feb75b13c08e2b057f25ac10e727e1ada81fd2b831628041ebb,
    curve,
)
# fmt: off
naf_pairs = [
    [
        1,
        [1]
    ],
    [
        3,
        [-1,0,1]
    ],
    [
        4,
        [0, 0, 1]
    ],
    [
        0b1111111111111111111111111111111111111111111111,
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1]
    ],
    [
        scalar,
        [-1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0,
        -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0,
        0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, -1,
        0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, -1, 0, -1,
        0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0,
        -1, 0, 0, 0, 0, 0, 1]
    ]
]
# fmt: on

print("[*] Teste NAF")
err_pairs = []
try:
    for binary, solution in naf_pairs:
        naf = ecc.calc_naf_representation(binary)
        if naf != solution:
            err_pairs.append((binary, solution, naf))
except Exception:
    print("[!] Crash")
    traceback.print_exc()
    print("[!] Ende")
    exit(-1)
if err_pairs:
    print("[-] NAF-Berechnung fehlerhaft")
    for binary, solution, naf in err_pairs:
        print(f"NAF for {binary} should be {solution} but was {naf}")
    print("[!] Ende")
    exit(-1)

print("[*] Teste Double-and-Add")
try:
    p_test = ecc.double_and_add(point, scalar)
    if p_test != compare_val:
        print("[-] Fehlerhaftes Ergebnis in Standard double-and-add")
        exit(-1)
    p_test = ecc.naf_double_and_add(point, scalar)
    if p_test != compare_val:
        print("[-] Fehlerhaftes Ergebnis in NAF double-and-add")
        exit(-1)
except Exception:
    print("[!] Crash")
    traceback.print_exc()
    print("[!] Ende")
    exit(-1)

random.seed(1)
for i in range(0, 100):
    factor = random.randint(0, 1<<255)
    if ecc.naf_double_and_add(point, factor) != ecc.double_and_add(point, factor):
        print("[-] Standard und NAF double-and-add verhalten sich unterschiedlich")
        print("[!] Ende")
        exit(-1)

print("[+] PASS")
