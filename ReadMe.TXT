�koly:
1) �prava videa - prokl�d�n�, rozlo�en� na sn�mky? (zrychlen� algoritmu?), filtrace?
2) N�vrh Gui
3) Vytvoren� modelu pozad� - Model prostred� poc�tat ve fin�le pomoc� medi�nu - robustnej�� ne� metoda prumerov�n�
5) Detekce j�zdn�ch pruhu 
6) Detekce pohybuj�c�ch se objektu - hranov� detektory, metoda rozd�lu modelu a sn�mku
7) Klasifikace automobilu - 
8) Optimalizace
9) Poc�t�n� ovecek
10) Dokumentace

Zdroje:
http://www.vutbr.cz/www_base/zav_prace_soubor_verejne.php?file_id=26890
http://www.mathworks.com/help/images/functionlist.html
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.84.9651&rep=rep1&type=pdf
http://www.inf.ed.ac.uk/teaching/courses/av/index.html

Probl�my:
- detekce pruhu - mu�e je definovat u�ivatel? Jak re�it, kdy� v pruhu projede jen nekolik m�lo aut -- myslim, ze je vyresene :)
- m� detekce aut prob�hat v re�ln�m case? -- Pochybuju, ze se to vubec da v matlabu implementovat

- Tahle posledni metoda je strasne vykonove narocna, co ji delat treba jednou za sec a pak ty stredy vyuzit pro shlukovou analyzu na jednoduchou derivaci snimku?
	- zkusit se to muze, ale nejdriv bych zkusil snizit rozliseni videa
	- znacne jsem zrychlil detekci, ale rychlejsi to uz asi nebude, chtelo by to CUDA nebo obdobu.
- v repu je stara verze get_background - script run nejde tedy spustit

- filtrovat fiktivni objekty


Probl�mov� �seky (snimky):

Video 00013.avi
930 - auto zmeni svuj centroid  a je zapocitano x2
1728 - auto zmeni svuj centroid  a je zapocitano x2, pak hodne aut
2500 - kamion pocitan x2