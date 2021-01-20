## How to update the website

### English version

- `https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/`以下のファイルを（`git clone`などで）ダウンロードしておく

- 以下、今週の更新日付を`2021XXXX`、先週の更新日付を`2021YYYY`とする（ファイル、フォルダ名に使用）

- 以下の説明にあるファイルおよびフォルダ名はダウンロードしたフォルダからの相対パス（たとえば、`/_archives/2021XXXX/`は、`(ダウンロードしたフォルダ名)/Covid19OutputJapan.github.io/_archives/2021XXXX/`）

0. `/_archives/2021XXXX/Main_Japan.m`をMATLABで実行して、Figure(.pngファイル)を`/image/2021XXXX`以下に保存

  - コードやデータはすべて`/_archives/2021XXXX`に保存

  - `/_archives/2021XXXX/Figure_JP.m`をMATLABで実行して、Figure(.pngファイル)を`/image/2021XXXX`以下に保存（日本語版サイトに使用）

1. `./2021XXXX.md`を`/index.md`からコピーして作成し、以下のように変更

  - 6行目の
  ```
  permalink: index.html
  ```
  を
  ```
  permalink: 2021XXXX.html
  ```
  に変更

  - （たとえば、`YYYY=0113`のとき）10行目の
  ```
  ## Updated weekly (Last update on January 13, 2021)
  ```
  を
  ```
  ## Updated on January 20, 2021
  ```
  に変更

2. `/index.md`を以下のように変更

  - （たとえば、`XXXX=0120`のとき）10行目を以下のように変更
  ```
  ## Updated weekly (Last update on January 20, 2021)
  ```

  - 図へのリンクをそれぞれ変更
    - 17行目
    ```
    |![Projection](./images/2021XXXX/VariablesProjection.png)|
    ```
    - 24行目
    ```
    |![TradeoffUB](./images/2021XXXX/BaselineTradeoffUB.png)<br>![Tradeoff](./images/2021XXXX/LaggedTradeoff.png)|
    ```
    - 33行目
    ```
    |![ForecastErrorsD](./images/2021XXXX/ForecastErrorsD.png)|
    ```
    - 38行目
    ```
    |![ForecastErrorsN](./images/2021XXXX/ForecastErrorsN.png)|
    ```
  - 表の数値（45-46,53-54行目）をそれぞれ変更(`Main_Japan.m`の結果を書き留めておく)

3. `/_data/sidebars/home_sidebar.yml`の23行目以降を次のように変更

  - `Latest`と`Last week`のリンクを変更し、2週前のリンクを加える
  
  - たとえば、`XXXX=0120`,`YYYY=0113`,`ZZZZ=0106`のとき
  ```
  - title: Nationwide
    output: web, pdf
    folderitems:
    - title: Latest
      url: /2021XXXX.html
      output: web, pdf
      type: homepage
    - title: Last week
      url: /2021YYYY.html
      output: web, pdf
    - title: January 6, 2021
      url: /2021ZZZZ.html
      output: web, pdf
  ```

4. 以下のgitコマンドを実行
```
git add -A
git commit -m "update on 2021XXXX"
git push
```

### Japanese version

- `https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/JP/`以下のファイルを（`git clone`などで）ダウンロードしておく
  - English versionと同じ構造になっています
  
0. 英語版サイトにある画像ファイルフォルダ`/Covid19OutputJapan.github.io/image/2021XXXX`を、日本語版サイト`/Covid19OutputJapan.github.io/JP/image/2021XXXX`に移動

1-4. 英語版サイトと同様の作業を行う
  - English versionのフォルダと混同しないように注意