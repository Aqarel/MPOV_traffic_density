Úkoly:
1) úprava videa - prokládání, rozložení na snímky? (zrychlení algoritmu?), filtrace?
2) Návrh Gui
3) Vytvorení modelu pozadí - Model prostredí pocítat ve finále pomocí mediánu - robustnejší než metoda prumerování
5) Detekce jízdních pruhu 
6) Detekce pohybujících se objektu - hranové detektory, metoda rozdílu modelu a snímku
7) Klasifikace automobilu - 
8) Optimalizace
9) Pocítání ovecek
10) Dokumentace

Zdroje:
http://www.vutbr.cz/www_base/zav_prace_soubor_verejne.php?file_id=26890

Problémy:
- detekce pruhu - muže je definovat uživatel? Jak rešit, když v pruhu projede jen nekolik málo aut
- má detekce aut probíhat v reálném case?
- okolo 216 snimku algritmus zkolabuje, pravdepodobne mrak, nebo neco to rozhodi. Treba si vic pohrat s prahama. Monza ten median pomuze.
- nebylo by lepsi delat detekci cistym rozdilem jednotlivych snimku a pak hodit prah?
- Tahle posledni metoda je strasne vykonove narocna, co ji delat treba jednou za sec a pak ty stredy vyuzit pro shlukovou analyzu na jednoduchou derivaci snimku?