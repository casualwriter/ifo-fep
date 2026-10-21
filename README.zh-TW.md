# ifo-fep：Markov blanket 的宣稱取決於觀測方式

[English](README.md) | **繁體中文**

CK Hung 與 Emmy Team（Amoy Technology）的兩篇相連的工作論文，以及程式碼和檢查。1.0 版，2026 年 10 月。

| 論文 | 內容 |
|---|---|
| [介面是與誰相接？Markov blanket 的宣稱取決於觀測方式](papers/fep-blanket/20261020-fep-blanket-v1.0-zh-tw.md)（[PDF](papers/fep-blanket/20261020-fep-blanket-v1.0-zh-tw.pdf)；[English](papers/fep-blanket/20261020-fep-blanket-v1.0.md)） | 主論文：觀點、主要結果，以及尚未證明的部分。 |
| [線性擴散中的有限時間窗路徑屏蔽：一個精確條件與殘餘漏量公式](papers/path-blanket/20261020-path-blanket-v1.0-zh-tw.md)（[PDF](papers/path-blanket/20261020-path-blanket-v1.0-zh-tw.pdf)；[English](papers/path-blanket/20261020-path-blanket-v1.0.md)） | 技術稿：證明、數值檢查，以及失敗預測與更正的揭露（附錄 A）。 |

## 一段話

自由能原理把 Markov blanket 說成內外之間的介面。我們在這個想法上再問一件事：是與誰相接的介面？對線性隨機系統，我們給出在不同觀測方式下 blanket 能否屏蔽內外的精確條件、把「屏蔽」和「編碼」分開、一個單檔工具，以及一個仍然無法判定的真實資料案例。這個觀點來自我們正在發展的更大框架（IFO）；所有結果不依賴它。我們希望專家能檢驗、修正，並帶它走得更遠。

## 程式碼

- `code/tools/blanket_report.py`：單一檔案，報告一個提議的 blanket 在不同觀察者眼中的樣子。不帶參數執行即為示範。需要 NumPy 和 SciPy。
- `code/track-a/`：每項檢查一支腳本（M1–M17d 及更早）。每支腳本在開頭寫明預測（或「事後」），印出結果，並寫出一個 JSON 檔（已附上）。技術稿的「可重現性」表對應腳本與章節。
- `code/tools/md2pdf.py`：把 Markdown 論文轉成 PDF（Chrome 與 KaTeX）。`make_release.py` 產生乾淨的發表版。`zenodo_upload.py` 上傳 Zenodo 草稿。

部分腳本開頭提到內部規劃筆記（`plan/...`、`inbox/...`），不在這個倉庫內；每支腳本開頭的預測與結果本身是完整的。

未附上：fMRI 資料（Gilson 等人 2018 的公開複本，github.com/mb-BCA/notebooks_review2019）、pyMOU（github.com/mb-BCA/pyMOU，commit f0626be）、pyDMBD（github.com/bayesianempirimancer/pyDMBD，commit e058e9c）。只有 fMRI 與 DMBD 的檢查需要它們。

## 授權與引用

論文：CC BY 4.0。程式碼：MIT。見 `LICENSE` 與 `LICENSE-papers.md`。引用方式見 `CITATION.cff`。1.0 版的 DOI：[10.5281/zenodo.23524009](https://doi.org/10.5281/zenodo.23524009)。所有版本的 DOI（永遠指向最新版）：[10.5281/zenodo.23524008](https://doi.org/10.5281/zenodo.23524008)。紀錄：https://zenodo.org/records/23524009

## 回饋

歡迎更正、反例，以及前人工作的線索：請開 issue。
