# pgx_recommender (v0.3)

CYP2C19 (+ opsiyonel CYP2D6) diplotipinden, CPIC rehberlerine göre 12 antidepresan (8 TCA/SSRI + 4 yalnızca-CYP2D6 ilaç) için sıralama üretir.
Bu sürümde **ML yok**; her adım CPIC tablolarından deterministik lookup.

## Kullanım

```bash
pip install pandas pytest
python -m pgxrec.build_data                      # ham xlsx -> data/processed/*.csv (bir kez)
python -m pytest -q                              # 53 test
python -m pgxrec.cli --cyp2c19 "*2/*17"
python -m pgxrec.cli --cyp2c19 "*1/*1" --cyp2d6 "*4/*4"
```

```python
from pgxrec.recommend import recommend, format_result
print(format_result(recommend("*1/*17", "*1/*4")))
```

## Akış

```
diplotip -> fenotip (CPIC tablosu) -> ilaç başına CPIC önerisi -> skor (0-3) -> sıralama
```

- **CYP2C19:** allel -> fonksiyon (Allele Functionality Table), fonksiyon çifti -> fenotip (CPIC fenotip tablosu).
  Tabloda olmayan allel "Uncertain function" sayılır ve işaretlenir (tahmin yürütülmez).
- **CYP2D6:** diplotip, CPIC'in hazır Diplotype-Phenotype tablosundan doğrudan okunur (kopya sayısı gösterimi dahil, örn. `*1/*2x2`).
  Allel sırası önemsizdir. Tabloda olmayan diplotip "Indeterminate".
- **Öneri tabloları:** (CYP2C19 fenotipi, CYP2D6 fenotipi) çiftiyle anahtarlanır; aktivite skoru aralıklarıyla uğraşılmaz.
  CYP2D6 verilmezse CPIC'in `CYP2D6 No Result` satırları kullanılır. Sertralin için `CYP2B6 No Result` satırları (CYP2B6 henüz desteklenmiyor).
- **Yalnızca CYP2D6'ya bağlı ilaçlar** (nortriptyline, desipramine, paroxetine, fluvoxamine): CPIC tablolarında "No Result" satırı yok.
  CYP2D6 verilmezse bu ilaçlar **sıralanmaz**, çıktıda "needs CYP2D6" olarak listelenir. CYP2C19'dan bağımsızdırlar
  (CYP2C19 Indeterminate olsa bile CYP2D6 verilirse sıralanırlar).
- **Aktivite skoruna duyarlı satır:** nortriptyline CYP2D6 IM'de öneri gücü aktivite skoruna göre değişir (0.25–0.5 Moderate, 0.75–1.0 Optional);
  build bu tür çelişkileri otomatik bulup aktivite skorunu anahtara ekler (başka ilaçta böyle bir durum yok).

## Skor (ilaç başına tek skor, düşük = iyi)

| Skor | Anlamı |
|---|---|
| 0 | standart doz |
| 1 | doz ayarı |
| 2 | avoid / alternatif (Optional güç) |
| 3 | avoid / alternatif (Moderate veya Strong güç) = **kritik** |
| - | CPIC "No recommendation": sıralanmaz, ayrıca listelenir |

Eşitlikler eşitlik olarak raporlanır. Öneri metni -> eylem eşlemesi `config/action_map.csv` içinde ve **elle gözden geçirilmiştir**;
rehber güncellenip yeni bir metin gelirse `build_data` hata verir (sessizce yanlış sınıflandırma olmaz).

## Bilinen limitasyonlar (tezde yazılmalı)

1. Skor etkililiği değil, CPIC'e göre **genetik uyumsuzluk riskini** ölçer.
2. "Kritik" sınırı CPIC'in öneri gücü etiketine (Strong/Moderate/Optional) bağlı. Bu etiket kısmen kanıt miktarını yansıtır:
   çok çalışılmış bir ilaç (örn. amitriptilin) daha "kritik" görünebilir, az çalışılmış ilaç daha güvenli görünebilir.
   Örnek: CYP2D6 PM'de amitriptilin Strong, diğer dört TCA Optional. Öneri: eylem-only skorla duyarlılık analizi yap.
3. Sertralin PM satırı "dozu azalt **ya da** alternatif seç" der; doz ayarı (1) sayıldı (karar, `action_map.csv`'de değiştirilebilir).
4. CYP2D6'da kopya sayısı/hibrit varyantlar standart VCF'ten güvenilir çağrılamaz; bu paket CYP2D6 diplotipini **girdi olarak** alır.
5. CYP2D6 verisi yoksa TCA sonuçları eksik bilgiyle üretilir (`[incomplete: no CYP2D6]` bayrağı) ve 4 yalnızca-CYP2D6 ilaç sıralanmaz:
   yani aday havuzu girdiye göre 8 ya da 12 ilaç olur, farklı girdilerin sıralamaları birbiriyle doğrudan karşılaştırılmamalı.
6. Popülasyon/ırk etkisi, yaş, eşlik eden ilaçlar (phenoconversion) modellenmedi.
7. Citalopram ve escitalopram CPIC'te birebir aynı tablodur; her zaman aynı grupta çıkar.

## Girdi ve referans genom (karar)

CPIC allele definition tabloları pozisyonları **GRCh38 (NC_000010.11, GRCh38.p2)** ve RefSeqGene (NG_008384.3) koordinatlarında verir.
Bu yüzden birincil girdi **GRCh38 VCF** olacak (faz bilgisi varsa kullanılır; yoksa belirsizlik bayrağı). FASTA yolu, diziyi GRCh38'e hizalayıp VCF üretecek ön modül olacak.
Dikkat: tabloda `*38` referans-dizi haplotipidir (g.94842866 pozisyonunda A), `*1` ise orada G taşır (rs3758581, p.I331V).
Yani GRCh38'e göre çağrılan bir VCF'te `*1` taşıyıcıları varyant (A>G) görünür; "varyant yok" = `*38`. Allel eşleme bu farkı doğru ele almalı
(`*38` ve `*1` ikisi de normal fonksiyon, fenotip değişmez, ama eşleme mantığı bozulursa herkes yanlış allele düşer).

## Doğrulama (gerçek veri)

`python -m pgxrec.validate`: NYGC 1000 Genomes GRCh38 faz VCF'inden (CYP2C19 bölgesi, 88 örnek) VCF -> allel eşleme (`pgxrec/matcher.py`) -> diplotip -> fenotip -> sıralama,
GeT-RM konsensüs genotipleriyle karşılaştırılır (`data/validation/`). Sonuç: diplotip, fenotip ve nihai öneri 88/88 uyumlu.

Dürüstlük notu: ilk çalıştırmada diplotip uyumu 61/88 idi (`*2` yerine `*35`, fenotip 87/88). Neden, kuralın tablodaki belirsiz (IUPAC) hücreleri saymamasıydı;
kural düzeltildi (belirsiz hücreler de tanımlı pozisyon sayılır, eşitlikte en düşük fonksiyonlu allel seçilir). Bu düzeltme aynı veri üzerinde yapıldığı için
88/88 bağımsız bir test sayılmaz; yeni örneklerle yeniden doğrulanmalı. 13 örnek (101'den) VCF'te yoktu. Nadir allel örnekleri çok az (`*3`, `*4`, `*8`, `*9`...).

## Sıradaki adımlar

- FASTA -> VCF ön modülü (yalnızca CYP2C19)
- Novel varyant modülü (missense, 16 etiketli allel; baseline karşılaştırmalı) -> "Uncertain function" yerine tahmini fonksiyon
- GeT-RM örnekleriyle doğrulama

Veri kaynağı: CPIC / ClinPGx referans tabloları (`data/raw`). Kullanırken CPIC rehberlerine atıf yapın.
