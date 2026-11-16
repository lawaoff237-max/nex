#!/usr/bin/env python3
"""
mask_search ni patternlar ustida ketma-ket ishga tushiradi.

Topilganni aniqlash mask_search chiqish formatiga BOG'LIQ EMAS:
chiqishdagi har bir hex qiymat va WIF Python ichida hash160 bilan tekshiriladi.

QAYTA ISHGA TUSHIRISHNI DAVOM ETTIRISH (resume) YO'Q: har bir ishga
tushirishda RO'YXATDAGI HAMMA patternlar boshidan qayta tekshiriladi.
Hech qanday "tekshirilgan.txt" fayli o'qilmaydi/yozilmaydi.

ASLIY PATTERNS FAYLI HECH QACHON JISMONAN BO'LINMAYDI: --patterns
doim XUDDI O'SHA bitta to'liq faylga ishora qiladi. Bo'lish faqat
MANTIQIY, --part / --parts orqali, xotirada amalga oshadi.

Ishlatish (Kaggle / parallel uchun argparse bilan):
  python3 tekshir_kaggle.py --patterns patterns.txt --part 1 --parts 100 --device 0
  python3 tekshir_kaggle.py --patterns patterns.txt --part 2 --parts 100 --device 1
  python3 tekshir_kaggle.py --patterns patterns.txt           -> --part berilmasa, HAMMASI
  python3 tekshir_kaggle.py --test     -> faqat 1-pattern, to'liq chiqish bilan (diagnostika)

--part X / --parts N : fayl X-1/N dan X/N gacha bo'lgan (dublikatsiz,
tartib saqlangan) qismga mantiqiy bo'linadi — xuddi split_patterns.py
--mode block bergan natija kabi, lekin diskka hech narsa yozilmasdan.
X qiymati 1 dan N gacha bo'lishi kerak.

Har bir --part uchun topilgan fayli avtomatik hosil qilinadi
(masalan --part 7 --parts 100 -> topilgan_part007.txt), shuning uchun
bir nechta qismni parallel ishga tushirish bir-biriga aralashmaydi.
"""
import argparse
import hashlib
import os
import re
import subprocess
import sys
from datetime import datetime

# ============================================================
# SOZLAMALAR (CLI orqali o'zgartiriladi, pastdagi parse_args() ga qarang)
# ============================================================
HASH160 = "f6f5431d25bbf7b12e8add9af5e3475c44a0a5b8"
PATTERNS_FILE = "patterns.txt"
TOPILGAN_FILE = "topilgan.txt"
MASK_SEARCH = "./mask_search"
DEVICE = None                # None = mask_search o'zining defaultini ishlatadi (0)
PART = None                  # None = bo'linmaydi, hammasi ishga tushadi
PARTS = 100
PATTERN_LENGTH = 71          # None qilsangiz uzunlik tekshirilmaydi
TIMEOUT_SEC = 86400


def parse_args():
    p = argparse.ArgumentParser(description="mask_search uchun ketma-ket ishga tushirish/tekshirish.")
    p.add_argument("--patterns", default=PATTERNS_FILE, help="to'liq patternlar fayli (hech qachon bo'linmaydi)")
    p.add_argument("--hash160", default=HASH160, help="qidirilayotgan hash160 (hex40)")
    p.add_argument("--mask-search", default=MASK_SEARCH, help="mask_search bajariladigan fayl yo'li")
    p.add_argument("--device", type=int, default=None, help="CUDA device raqami (--device argumentini mask_search ga uzatadi)")
    p.add_argument("--part", type=int, default=None,
                   help="mantiqiy qism raqami, 1..--parts (berilmasa: faylning hammasi ishga tushadi)")
    p.add_argument("--parts", type=int, default=PARTS,
                   help=f"jami mantiqiy qismlar soni (default: {PARTS})")
    p.add_argument("--suffix", default=None,
                   help="topilgan fayli uchun qo'shimcha (default: --part dan yoki --patterns nomidan avtomatik olinadi)")
    p.add_argument("--test", action="store_true", help="faqat 1-pattern, to'liq chiqish bilan (diagnostika)")
    return p.parse_args()

# ============================================================
# secp256k1 + hash160 (tekshiruv uchun)
# ============================================================
P = 2**256 - 2**32 - 977
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)

def _add(a, b):
    if a is None: return b
    if b is None: return a
    if a[0] == b[0] and (a[1] + b[1]) % P == 0: return None
    if a == b:
        l = 3 * a[0] * a[0] * pow(2 * a[1], -1, P) % P
    else:
        l = (b[1] - a[1]) * pow(b[0] - a[0], -1, P) % P
    x = (l * l - a[0] - b[0]) % P
    return (x, (l * (a[0] - x) - a[1]) % P)

def _mul(k, pt=G):
    r = None
    while k:
        if k & 1: r = _add(r, pt)
        pt = _add(pt, pt)
        k >>= 1
    return r

def _rol(x, n): return ((x << n) | (x >> (32 - n))) & 0xFFFFFFFF

def _ripemd160(data):
    try:
        h = hashlib.new("ripemd160"); h.update(data); return h.digest()
    except Exception:
        pass
    r1 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,7,4,13,1,10,6,15,3,12,0,9,5,2,14,11,8,
          3,10,14,4,9,15,8,1,2,7,0,6,13,11,5,12,1,9,11,10,0,8,12,4,13,3,7,15,14,5,6,2,
          4,0,5,9,7,12,2,10,14,1,3,8,11,6,15,13]
    r2 = [5,14,7,0,9,2,11,4,13,6,15,8,1,10,3,12,6,11,3,7,0,13,5,10,14,15,8,12,4,9,1,2,
          15,5,1,3,7,14,6,9,11,8,12,2,10,0,4,13,8,6,4,1,3,11,15,0,5,12,2,13,9,7,10,14,
          12,15,10,4,1,5,8,7,6,2,13,14,0,3,9,11]
    s1 = [11,14,15,12,5,8,7,9,11,13,14,15,6,7,9,8,7,6,8,13,11,9,7,15,7,12,15,9,11,7,13,12,
          11,13,6,7,14,9,13,15,14,8,13,6,5,12,7,5,11,12,14,15,14,15,9,8,9,14,5,6,8,6,5,12,
          9,15,5,11,6,8,13,12,5,12,13,14,11,8,5,6]
    s2 = [8,9,9,11,13,15,15,5,7,7,8,11,14,14,12,6,9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
          9,7,15,11,8,6,6,14,12,13,5,14,13,13,7,5,15,5,8,11,14,14,6,14,6,9,12,9,12,5,15,8,
          8,5,12,9,12,5,14,6,8,13,6,5,15,13,11,11]
    K1 = [0, 0x5A827999, 0x6ED9EBA1, 0x8F1BBCDC, 0xA953FD4E]
    K2 = [0x50A28BE6, 0x5C4DD124, 0x6D703EF3, 0x7A6D76E9, 0]
    def f(j, x, y, z):
        if j < 16: return x ^ y ^ z
        if j < 32: return (x & y) | (~x & z)
        if j < 48: return (x | ~y) ^ z
        if j < 64: return (x & z) | (y & ~z)
        return x ^ (y | ~z)
    m = bytearray(data); l = len(data) * 8
    m.append(0x80)
    while len(m) % 64 != 56: m.append(0)
    m += l.to_bytes(8, "little")
    h = [0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0]
    for off in range(0, len(m), 64):
        X = [int.from_bytes(m[off+4*i:off+4*i+4], "little") for i in range(16)]
        a, b, c, d, e = h; a2, b2, c2, d2, e2 = h
        for j in range(80):
            t = (_rol((a + f(j, b, c, d) + X[r1[j]] + K1[j//16]) & 0xFFFFFFFF, s1[j]) + e) & 0xFFFFFFFF
            a, e, d, c, b = e, d, _rol(c, 10), b, t
            t = (_rol((a2 + f(79-j, b2, c2, d2) + X[r2[j]] + K2[j//16]) & 0xFFFFFFFF, s2[j]) + e2) & 0xFFFFFFFF
            a2, e2, d2, c2, b2 = e2, d2, _rol(c2, 10), b2, t
        t = (h[1] + c + d2) & 0xFFFFFFFF
        h[1] = (h[2] + d + e2) & 0xFFFFFFFF
        h[2] = (h[3] + e + a2) & 0xFFFFFFFF
        h[3] = (h[4] + a + b2) & 0xFFFFFFFF
        h[4] = (h[0] + b + c2) & 0xFFFFFFFF
        h[0] = t
    return b"".join(x.to_bytes(4, "little") for x in h)

def hash160s(k):
    """Kalit uchun (siqilgan, siqilmagan) hash160 hex."""
    x, y = _mul(k)
    comp = bytes([2 + (y & 1)]) + x.to_bytes(32, "big")
    unc = b"\x04" + x.to_bytes(32, "big") + y.to_bytes(32, "big")
    return tuple(_ripemd160(hashlib.sha256(p).digest()).hex() for p in (comp, unc))


# ============================================================
# Base58Check va WIF kodlash / ochish
# ============================================================
BASE58_ALPHABET = "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz"

def _base58_encode(b: bytes) -> str:
    n = int.from_bytes(b, "big")
    res = ""
    while n > 0:
        n, r = divmod(n, 58)
        res = BASE58_ALPHABET[r] + res
    n_leading_zeros = len(b) - len(b.lstrip(b"\x00"))
    return "1" * n_leading_zeros + res

def _base58_decode(s: str):
    if not s or any(c not in BASE58_ALPHABET for c in s):
        return None
    n = 0
    for c in s:
        n = n * 58 + BASE58_ALPHABET.index(c)
    pad = 0
    for c in s:
        if c == "1":
            pad += 1
        else:
            break
    full = b"" if n == 0 else n.to_bytes((n.bit_length() + 7) // 8, "big")
    return b"\x00" * pad + full

def privkey_to_wif(k_int, compressed=True, mainnet=True):
    """Private key (int) -> WIF satr."""
    prefix = b"\x80" if mainnet else b"\xef"
    payload = prefix + k_int.to_bytes(32, "big") + (b"\x01" if compressed else b"")
    checksum = hashlib.sha256(hashlib.sha256(payload).digest()).digest()[:4]
    return _base58_encode(payload + checksum)

def wif_to_privkey(s: str):
    """WIF -> (kalit_int, 'compressed'|'uncompressed') yoki None."""
    raw = _base58_decode(s.strip())
    if raw is None or len(raw) < 5:
        return None
    payload, checksum = raw[:-4], raw[-4:]
    calc = hashlib.sha256(hashlib.sha256(payload).digest()).digest()[:4]
    if calc != checksum or not payload or payload[0] != 0x80:
        return None
    body = payload[1:]
    if len(body) == 32:
        return int.from_bytes(body, "big"), "uncompressed"
    if len(body) == 33 and body[-1] == 0x01:
        return int.from_bytes(body[:-1], "big"), "compressed"
    return None


# ============================================================
# Yordamchi funksiyalar
# ============================================================
def load_patterns():
    """To'liq faylni o'qiydi (hech qachon jismonan bo'linmaydi).
    PART berilgan bo'lsa, faqat o'sha mantiqiy bo'lakni qaytaradi:
    dublikatlar olib tashlanadi (tartib saqlanadi), so'ng fayl PARTS ga
    teng bo'laklarga bo'linadi (split_patterns.py --mode block bilan bir xil
    qoida) va faqat PART-chi bo'lak qaytariladi — diskka hech narsa yozilmaydi."""
    if not os.path.exists(PATTERNS_FILE):
        print(f"XATO: {PATTERNS_FILE} topilmadi!")
        sys.exit(1)
    with open(PATTERNS_FILE, encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

    if PART is None:
        return lines

    seen = set()
    unique = []
    for l in lines:
        if l not in seen:
            seen.add(l)
            unique.append(l)

    if PARTS < 1:
        print(f"XATO: --parts musbat bo'lishi kerak (olingan: {PARTS})")
        sys.exit(1)
    if not (1 <= PART <= PARTS):
        print(f"XATO: --part {PARTS} ta qismdan biri bo'lishi kerak, ya'ni 1..{PARTS} (olingan: {PART})")
        sys.exit(1)

    size = (len(unique) + PARTS - 1) // PARTS
    start = (PART - 1) * size
    end = start + size
    return unique[start:end]


def _append_line(path, text):
    with open(path, "a", encoding="utf-8") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())


def save_found(pattern, key_hex, matched_format, wif_compressed, wif_uncompressed):
    wif_mos = wif_compressed if matched_format == "compressed" else wif_uncompressed
    _append_line(
        TOPILGAN_FILE,
        f"{datetime.now()} | pattern={pattern} | key_hex={key_hex} | "
        f"mos_format={matched_format} | WIF_mos={wif_mos} | "
        f"WIF_compressed={wif_compressed} | WIF_uncompressed={wif_uncompressed}\n",
    )


def _chunks_to_text(chunks):
    parts = []
    for chunk in chunks:
        if not chunk:
            continue
        if isinstance(chunk, bytes):
            parts.append(chunk.decode(errors="ignore"))
        else:
            parts.append(chunk)
    return "".join(parts)


def run_mask_search(pattern):
    """Qaytaradi: (cmd, returncode, chiqish, davomiylik, holat).
    Hech qanday log fayl yaratmaydi — faqat xotiradagi output qaytaradi."""
    cmd = [MASK_SEARCH, "--mask", pattern, "--hash160", HASH160]
    if DEVICE is not None:
        cmd += ["--device", str(DEVICE)]
    start = datetime.now()
    rc, status = None, "ok"
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_SEC)
        rc = r.returncode
        output = (r.stdout or "") + (r.stderr or "")
    except subprocess.TimeoutExpired as e:
        output = "TIMEOUT\n" + _chunks_to_text((e.stdout, e.stderr))
        status = "timeout"
    except FileNotFoundError:
        print(f"XATO: {MASK_SEARCH} topilmadi (yo'l yoki chmod +x tekshiring)!")
        sys.exit(1)
    except Exception as e:
        output = f"ERROR: {e}"
        status = "error"
    return cmd, rc, output, (datetime.now() - start).total_seconds(), status


def _check_key(k, seen, target):
    if k in seen or not (1 <= k < N):
        return None
    seen.add(k)
    comp_hash, unc_hash = hash160s(k)
    if target == comp_hash:
        return format(k, "064x"), "compressed"
    if target == unc_hash:
        return format(k, "064x"), "uncompressed"
    return None


def verify_keys(output):
    """Chiqishdagi hex va WIF qiymatlarni private key deb olib, hash160 ni tekshiradi.
    Mos kelsa (kalit_hex, qaysi_format_mos_keldi) qaytaradi, aks holda (None, None)."""
    target = HASH160.lower()
    seen = set()

    for tok in re.findall(
        r"(?<![1-9A-HJ-NP-Za-km-z])[5KL][1-9A-HJ-NP-Za-km-z]{50,51}(?![1-9A-HJ-NP-Za-km-z])",
        output,
    ):
        parsed = wif_to_privkey(tok)
        if not parsed:
            continue
        hit = _check_key(parsed[0], seen, target)
        if hit:
            return hit

    for tok in re.findall(r"[0-9a-fA-F]{4,64}", output):
        try:
            k = int(tok, 16)
        except ValueError:
            continue
        hit = _check_key(k, seen, target)
        if hit:
            return hit
    return None, None


def claims_real_found(output):
    """Chiqish haqiqatan ham 'topildi' deyaptimi?
    not-found / notfound / found: 0 kabi inkorlar hisobga olinmaydi."""
    low = output.lower()
    low = re.sub(r"not[\s_\-]*found", " ", low)
    low = re.sub(r"\bnothing\s+found\b", " ", low)
    low = re.sub(r"\bno\s+keys?\s+found\b", " ", low)
    low = re.sub(r"\b0\s+keys?\s+found\b", " ", low)
    low = re.sub(r"\b0\s+found\b", " ", low)
    low = re.sub(r"\bfound\s*[:=]?\s*0\b", " ", low)
    low = re.sub(r"\bzero\s+found\b", " ", low)
    return re.search(r"\bfound\b", low) is not None


def decide(status, rc, output):
    """Kalit tasdiqlanmagan chiqishni tasniflaydi (faqat log/statistika uchun)."""
    if status == "timeout":
        return "timeout"
    if status == "error":
        return "error"
    # rc=0 = FOUND. Kalit tasdiqlanmasa, tekshirilgan deb yozib bo'lmaydi.
    if rc == 0:
        return "rc_found_unverified"
    if rc != 1:
        return "bad_rc"
    if not output.strip():
        return "empty"
    if claims_real_found(output):
        return "text_found_unverified"
    return "notfound"


# ============================================================
# Asosiy qism
# ============================================================
def main():
    global HASH160, PATTERNS_FILE, TOPILGAN_FILE, MASK_SEARCH, DEVICE, PART, PARTS

    args = parse_args()
    HASH160 = args.hash160
    PATTERNS_FILE = args.patterns
    MASK_SEARCH = args.mask_search
    DEVICE = args.device
    PART = args.part
    PARTS = args.parts

    # topilgan faylini --part dan (yoki --patterns nomidan) hosil qilish,
    # shunda bir nechta qism (--part 1, --part 2, ...) bir-biriga aralashmaydi.
    suffix = args.suffix
    if suffix is None:
        if PART is not None:
            width = len(str(PARTS))
            suffix = f"_part{PART:0{width}d}"
        else:
            base = os.path.basename(PATTERNS_FILE)
            stem, _ext = os.path.splitext(base)
            suffix = stem.replace("patterns", "", 1) if stem.startswith("patterns") else f"_{stem}"
    TOPILGAN_FILE = f"topilgan{suffix}.txt"

    test_mode = args.test
    patterns = load_patterns()
    if test_mode:
        patterns = patterns[:1]

    print("=" * 60)
    print(f"Boshlangan: {datetime.now()}")
    print(f"Hash160: {HASH160}")
    print(f"Patterns fayli: {PATTERNS_FILE} (jismoniy bo'linmagan)")
    if PART is not None:
        print(f"Mantiqiy qism: {PART}/{PARTS}")
    print(f"Device: {DEVICE if DEVICE is not None else '(default)'}")
    print(f"Topilgan fayli: {TOPILGAN_FILE}")
    print(f"Jami patternlar (shu qismda): {len(patterns)} (qayta ishga tushirishda HAMMASI qayta tekshiriladi, resume yo'q)")
    if test_mode:
        print("TEST REJIMI: faqat 1-pattern, to'liq chiqish ko'rsatiladi")
    print("=" * 60)

    total = len(patterns)
    found_count = problems = 0

    try:
        for i, pattern in enumerate(patterns, 1):
            if PATTERN_LENGTH and len(pattern) != PATTERN_LENGTH:
                print(f"[{i}/{total}] ⚠️  Uzunlik xato: {len(pattern)} ({PATTERN_LENGTH} kutilgan) — o'tkazildi")
                problems += 1
                continue

            print(f"\n[{i}/{total}] Pattern: {pattern}")
            cmd, rc, output, duration, status = run_mask_search(pattern)

            tail_lines = [l for l in output.strip().splitlines() if l.strip()]
            shown = tail_lines if test_mode else tail_lines[-3:]
            for l in shown:
                print(f"  | {l}")

            speeds = re.findall(r"([\d.]+)\s*Mkey/s", output)
            if speeds:
                print(f"  ⚡️ Tezlik: {speeds[-1]} Mkey/s")

            key, matched_format = verify_keys(output)
            if key:
                k_int = int(key, 16)
                wif_c = privkey_to_wif(k_int, compressed=True)
                wif_u = privkey_to_wif(k_int, compressed=False)
                wif_mos = wif_c if matched_format == "compressed" else wif_u
                print(f"  ✅ TOPILDI va TASDIQLANDI! Kalit (hex): {key}")
                print(f"  ✅ Mos kelgan format: {matched_format}")
                print(f"  🔑 WIF (mos, shu manzil): {wif_mos}")
                print(f"  🔑 WIF (compressed): {wif_c}")
                print(f"  🔑 WIF (uncompressed): {wif_u}")
                print(f"  ⏱️  Vaqt: {duration:.2f}s")
                save_found(pattern, key, matched_format, wif_c, wif_u)
                found_count += 1
                break

            outcome = decide(status, rc, output)
            if outcome == "timeout":
                print(f"  ⏰ TIMEOUT ({duration:.2f}s)")
                problems += 1
                continue
            if outcome == "error":
                print(f"  ⚠️  XATO ({duration:.2f}s)")
                problems += 1
                continue
            if outcome == "rc_found_unverified":
                print(f"  ⚠️  rc=0 (FOUND), lekin kalit hash160 ga mos kelmadi ({duration:.2f}s)")
                problems += 1
                continue
            if outcome == "bad_rc":
                print(f"  ⚠️  Kutilmagan rc={rc} ({duration:.2f}s)")
                problems += 1
                continue
            if outcome == "text_found_unverified":
                print(f"  ⚠️  Chiqishda 'found' bor, lekin kalit hash160 ga mos kelmadi!")
                problems += 1
                continue
            if outcome == "empty":
                print(f"  ⚠️  Bo'sh chiqish ({duration:.2f}s) — noaniq")
                problems += 1
                continue

            print(f"  ❌ NOTFOUND ({duration:.2f}s)")

    except KeyboardInterrupt:
        print("\n⛔️ Foydalanuvchi tomonidan to'xtatildi.")

    print("\n" + "=" * 60)
    print(f"Tugadi: {datetime.now()}")
    print(f"Jami: {total} | Topildi: {found_count} | Muammoli: {problems}")
    print("=" * 60)


if __name__ == "__main__":
    main()