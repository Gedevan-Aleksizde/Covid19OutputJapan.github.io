## How to update the website

Last updated: January 22, 2021, written by Takeki Sunakawa

### English version

#### 0. 必要なソフトウェアのインストール、その他の注意点

- (Windowsの場合)git for windowsを[ここ](https://gitforwindows.org/)からダウンロードする。インストールするとき、改行コードの自動変換指定は"checkout as-is, commit unix-style line endings"を選択。あとはすべてデフォルトでよい

- (Macの場合)ターミナルから`git --version`でインストールされているか確認

- 初めてgitを使う場合、

  ```
  git config --global user.name "FIRST_NAME LAST_NAME"
  git config --global user.email "MY_NAME@example.com"
  ```
  として名前、メールアドレスを設定

- `https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/`以下のファイルを

  ```
  git clone https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/
  ```
  などとしてダウンロードしておく

- マークダウンファイル(.mdファイル)を編集する際、文字コードは`utf-8, ja-JP`とすること。特にWindowsの場合に注意。

- 以下の説明で、今週の更新日付を`2021XXXX`、先週の更新日付を`2021YYYY`とする（ファイル、フォルダ名に使用）。また、ファイルおよびフォルダ名はダウンロードしたフォルダからの相対パス（たとえば、`/_archives/2021XXXX/`は、`(ダウンロードしたフォルダ名)/Covid19OutputJapan.github.io/_archives/2021XXXX/`）

#### 1. MATLABを実行してFigureを作成

- `/_archives/2021XXXX/Main_Japan.m`をMATLABで実行して、Figure(.pngファイル)を`/image/2021XXXX`以下に保存

- コードやデータはすべて`/_archives/2021XXXX`に保存

<!--  - `/_archives/2021XXXX/Figure_JP.m`をMATLABで実行して、Figure(.pngファイル)を`/image/2021XXXX`以下に保存（日本語版サイトに使用） -->

#### 2. 先週のトップページ(`index.md`)を`2021YYYY.md`としてコピー、編集して保存

- `/2021YYYY.md`を`/index.md`からコピーして作成し、以下のように変更

  - 6行目の
  ```
  permalink: index.html
  ```
  を
  ```
  permalink: 2021YYYY.html
  ```
  に変更

#### 3. 今週の`index.md`を作成

- `/index.md`を以下のように変更

  - （たとえば、`XXXX=0120`のとき）10行目?を以下のように変更
  ```
  ## Updated on January 20, 2021
  ```

  - 図へのリンクをそれぞれ変更。たとえば、
  ```
  |![Projection](./images/2021XXXX/VariablesProjection.png)|
  ```
  
  - 表の数値をそれぞれ変更(`Main_Japan.m`の結果を書き留めておく)。たとえば、
  ```
  | **New Cases** | 53,088   |  41,290  | <span style="color: black; ">11,798</span> |
  | **New Deaths** |   723  | 445  | <span style="color: black; ">278</span> |
  ```
  ```
  | **New Cases** |  83,138  |  129,454  | <span style="color: red; ">-46,315</span> |
  | **New Deaths** |   1,004  |    1,459 | <span style="color: red; ">-454</span> |
  ```
    - 予測誤差がマイナスの場合は、フォントの色を赤にする
    
#### 4. サイドバー(`/_data/sidebars/home_sidebar.yml`)の更新

- `/_data/sidebars/home_sidebar.yml`の23行目?以降を次のように変更

  - `Latest`と`Last week`のリンクを変更し、2週前のリンクを加える

  - たとえば、`XXXX=0120`, `YYYY=0113`, `ZZZZ=0106`のとき
  ```
  - title: Nationwide
    output: web, pdf
    folderitems:
    - title: Latest
      url: /20210120.html
      output: web, pdf
      type: homepage
    - title: Last week
      url: /20210113.html
      output: web, pdf
    - title: January 6, 2021
      url: /20210106.html
      output: web, pdf
  ```

#### 5. アップロード

- 以下のgitコマンドを実行
```
git add -A
git commit -m "update on 2021XXXX"
git push
```

- （特にWindowsの場合）`git add -A`で時間がかかる場合、`git add (更新したファイル)`とする

- 初回の`git push`で、githubのCovid19OutputJapanアカウントで認証手続きが必要な場合がある

### Japanese version

#### 0.
- `https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/JP/`以下のファイルを（`git clone`などで）ダウンロードしておく
  - English versionと同じ構造になっています

#### 1. 
- 英語版サイトにある画像ファイルフォルダ`/Covid19OutputJapan.github.io/image/2021XXXX`を、日本語版サイト`/Covid19OutputJapan.github.io/JP/image/2021XXXX`に移動またはコピー

#### 2-5. 
- 英語版サイトと同様の作業を行う
  - English versionのフォルダと混同しないように注意
