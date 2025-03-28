# MCP Simple PubMed (Improved Version)
[![smithery badge](https://smithery.ai/badge/mcp-simple-pubmed)](https://smithery.ai/server/mcp-simple-pubmed)

> This is a fork of [andybrandt/mcp-simple-pubmed](https://github.com/andybrandt/mcp-simple-pubmed) with enhanced functionality and improvements.

> これは[andybrandt/mcp-simple-pubmed](https://github.com/andybrandt/mcp-simple-pubmed)をフォークし、機能強化と改良を加えたバージョンです。

An MCP server that provides access to PubMed articles through the Entrez API.

<a href="https://glama.ai/mcp/servers/5wlfb8i6bj"><img width="380" height="200" src="https://glama.ai/mcp/servers/5wlfb8i6bj/badge" alt="mcp-simple-pubmed MCP server" /></a>

## Features | 機能

- Search PubMed database using keywords | キーワードによるPubMedデータベース検索
- Access article abstracts | 論文アブストラクトの取得
- Download full text when available (enhanced functionality) | 利用可能な場合のフルテキストダウンロード（機能強化）
- Option to convert XML to plain text | XMLからプレーンテキストへの変換オプション
- Improved handling of large papers | 大きな論文の処理機能の改善
- Enhanced search query processing | 検索クエリ処理の強化
- Additional metadata in search results | 検索結果に追加のメタデータを含む

---

## Enhancements Over Original Version | オリジナルバージョンからの改良点

### Full Text Retrieval Enhancements | フルテキスト取得機能の改善

- Integrated functionality of `fulltext_client.py` and `pubmed_fetch.py` to eliminate code duplication
- Improved robustness of PMC ID retrieval logic (supports multiple retrieval methods)
- Enhanced processing of large papers
- Added option to convert XML format to plain text
- Enriched fallback information when full text is not available

> `fulltext_client.py`と`pubmed_fetch.py`の機能を統合し、コードの重複を削除
> PMC ID取得ロジックの堅牢性を向上（複数の取得方法をサポート）
> 大きな論文のフルテキスト取得処理を改善
> XML形式からプレーンテキストへの変換オプションを追加
> フルテキストが利用できない場合の代替情報を充実

### Search Functionality Improvements | 検索機能の改善

- Enhanced search query processing (especially date range formatting)
- Added additional metadata to search results (MeSH terms, keywords, etc.)
- Strengthened search error handling
- Improved generation of paper resource URIs and URLs

> 検索クエリ処理の改善（特に日付範囲のフォーマット）
> 検索結果に追加メタデータを含める（MeSH用語、キーワードなど）
> 検索エラーハンドリングの強化
> 論文リソースURIとURLの生成を改善

### Server Functionality Enhancements | サーバー機能の改善

- Improved robustness of tool invocation processing
- Enhanced error messages and alternative information
- Added detailed logging

> ツール呼び出し処理の堅牢性を向上
> エラーメッセージと代替情報の充実
> 詳細なロギングを追加

---

## Usage Tips | 使用方法のヒント

### Full Text Format Options | フルテキスト形式オプション

You can now choose between XML format and text format when retrieving full text. XML format preserves document structure information, but text format is more human-readable.

> 全文取得時にXML形式とテキスト形式を選択できるようになりました。XML形式は文書構造情報を保持しますが、人間が読みやすいのはテキスト形式です。

### Improved Search Capabilities | 改善された検索機能

The following features have been improved in search queries:
- Proper handling of date ranges (e.g., `2020:2024[Date - Publication]`)
- Increased robustness of field-specific searches
- Support for retrieving up to 50 search results (default is 10)

> 検索クエリで以下の機能が改善されました：
> - 日付範囲の適切な処理（例：`2020:2024[Date - Publication]`）
> - フィールド特定検索の堅牢性向上
> - 50件までの検索結果取得に対応（デフォルトは10件）

### Fallback Information | フォールバック情報

When full text is not available, the following alternative information is provided:
- Link to PubMed page
- Link to publisher site via DOI
- Paper citation information and title

> 全文が利用できない場合、以下の代替情報が提供されます：
> - PubMedページへのリンク
> - DOI経由の出版社サイトへのリンク
> - 論文の引用情報とタイトル

---

## Original Features (From Base Project) | 元のプロジェクトの機能

Please note that the tool returns XML-ized version of full text. It is however more useful for AIs than a "human readable" text would have been as it gives them additional information about document's structure. At least, this is what Claude 3.5 Sonnet said he prefers. 

Please also note that inability of this tool and possibly other tools to deliver a paper's full text may not be due to the fact that it is not available. When testing this tool I came across a paper that did not have full text on PubMed and when Claude accessed the publication URL (which we did get through DOI) using fetch he did get a "forbidden" error. However, I was able to access the very same page using a regular browser. 

In other words if your AI assistant is not able to get the full text of a paper using this tool it is worth trying manually with a regular web browser.

Finally, this tool of course can't give you access to paywalled/paid papers. You may be able to read them through your library access or – as a last resort – through a certain site that strives to make publicly funded research freely available. 

## Installation | インストール方法

### Development Installation | 開発版インストール

If you want to use this enhanced version with all improvements, you can install it using either method below.

#### For Git Users | Gitユーザー向け

##### Clone the repository | リポジトリのクローン
```bash
git clone https://github.com/masa061580/mcp-simple-pubmed.git
cd mcp-simple-pubmed
```

##### Install in development mode | 開発モードでインストール
```bash
pip install biopython mcp
pip install -e .
```

#### For Beginners (No Git Required) | 初心者向け（Git不要）

こちらの方法は、Gitを使用せずにインストールする方法です。

This method allows you to install without using Git.

##### Step 0: Install Python (if not already installed) | Pythonのインストール（まだない場合）

If you don't have Python installed, follow these steps:

1. Download Python from the official website: https://www.python.org/downloads/
   - Click on "Download Python 3.x.x" (the latest version)

2. Run the installer:
   - **Important**: Check the box that says "Add Python to PATH" during installation
   - Click "Install Now"

3. Verify installation:
   - Open Command Prompt (Windows) or Terminal (Mac/Linux)
   - Type `python --version` and press Enter
   - You should see something like `Python 3.x.x`

> Pythonがまだインストールされていない場合は、以下の手順に従ってください：
>
> 1. Pythonを公式サイトからダウンロード：https://www.python.org/downloads/
>    - 「Download Python 3.x.x」（最新バージョン）をクリック
>
> 2. インストーラーを実行：
>    - **重要**: インストール中に「Add Python to PATH」のボックスにチェックを入れる
>    - 「Install Now」をクリック
>
> 3. インストールの確認：
>    - コマンドプロンプト（Windows）またはターミナル（Mac/Linux）を開く
>    - `python --version`と入力してEnterを押す
>    - `Python 3.x.x`のように表示されれば成功です

##### Step 1: Download ZIP file | ZIPファイルのダウンロード

1. Go to the GitHub repository page: https://github.com/masa061580/mcp-simple-pubmed
2. Click the green "Code" button
3. Select "Download ZIP" from the dropdown menu
4. Save the ZIP file to a location you can easily find (e.g., Desktop or Downloads folder)

> 1. GitHubリポジトリページにアクセス：https://github.com/masa061580/mcp-simple-pubmed
> 2. 緑色の「Code」ボタンをクリック
> 3. ドロップダウンメニューから「Download ZIP」を選択
> 4. ZIPファイルを簡単に見つけられる場所（例：デスクトップやダウンロードフォルダ）に保存

##### Step 2: Extract the ZIP file | ZIPファイルの解凍

1. Right-click on the downloaded ZIP file
2. Select "Extract All..." (Windows) or use your preferred extraction tool
3. Choose a destination folder (Remember this location!)
4. Click "Extract"

> 1. ダウンロードしたZIPファイルを右クリック
> 2. 「すべて展開...」（Windows）を選択、または好みの解凍ツールを使用
> 3. 解凍先フォルダを選択（この場所を覚えておいてください！）
> 4. 「展開」をクリック

##### Step 3: Install required packages | 必要なパッケージのインストール

1. Open Command Prompt (Windows) or Terminal (Mac/Linux)
   - Windows: Press `Windows key + R`, type `cmd` and press Enter
   - Mac: Open Finder > Applications > Utilities > Terminal

2. Install the required packages:
   ```bash
   pip install biopython mcp
   ```

> 1. コマンドプロンプト（Windows）またはターミナル（Mac/Linux）を開く
>    - Windows: `Windowsキー + R`を押し、`cmd`と入力してEnterを押す
>    - Mac: Finder > アプリケーション > ユーティリティ > ターミナルを開く
>
> 2. 必要なパッケージをインストール：
>    ```bash
>    pip install biopython mcp
>    ```

##### Step 4: Install the package in development mode | 開発モードでパッケージをインストール

1. In Command Prompt/Terminal, navigate to the extracted folder:
   ```bash
   cd path\to\extracted\mcp-simple-pubmed-master\mcp-simple-pubmed-master
   ```
   Replace `path\to\extracted\mcp-simple-pubmed` with the actual path where you extracted the ZIP file.
   
   Example (Windows):
   ```bash
   cd C:\Users\YourUsername\Downloads\mcp-simple-pubmed-master\mcp-simple-pubmed-master
   ```
   
   Example (Mac/Linux):
   ```bash
   cd /Users/YourUsername/Downloads/mcp-simple-pubmed-master/mcp-simple-pubmed-master
   ```

2. Install the package in development mode:
   ```bash
   pip install -e .
   ```

> 1. コマンドプロンプト/ターミナルで、解凍したフォルダに移動します：
>    ```bash
>    cd 解凍先パス\mcp-simple-pubmed
>    ```
>    `解凍先パス\mcp-simple-pubmed`を実際にZIPファイルを解凍した場所のパスに置き換えてください。
>    
>    例（Windows）：
>    ```bash
>    cd C:\Users\YourUsername\Downloads\mcp-simple-pubmed-master\mcp-simple-pubmed-master
>    ```
>    
>    例（Mac/Linux）：
>    ```bash
>    cd /Users/YourUsername/Downloads/mcp-simple-pubmed-master/mcp-simple-pubmed-master
>    ```
>
> 2. 開発モードでパッケージをインストール：
>    ```bash
>    pip install -e .
>    ```

This development mode installation allows you to use the code with all enhancements and make further modifications without reinstalling.

> 開発モードでインストールすると、すべての機能強化が含まれたコードを使用でき、さらに変更を加えても再インストールが不要になります。

### Installing via Smithery | Smitheryを使ったインストール
> **注意**: SmitheryとPipでインストールされるのは元のバージョンであり、この強化版ではありません。この改良版を使用するには、上記の開発版インストール方法を使用してください。

To install the original version of Simple PubMed for Claude Desktop automatically via [Smithery](https://smithery.ai/server/mcp-simple-pubmed):

```bash
npx -y @smithery/cli install mcp-simple-pubmed --client claude
```

### Manual Installation of Released Version | リリース版の手動インストール
```bash
pip install mcp-simple-pubmed
```

**Note**: The Smithery and pip installations will install the original version, not this enhanced fork. Use the development installation method above to use this improved version.



## Configuration | 設定

The server requires the following environment variables:

- `PUBMED_EMAIL`: Your email address (required by NCBI)
- `PUBMED_API_KEY`: Optional API key for higher rate limits 

The standard rate limit is 3 requests / second. No rate limiting was implemented, as it is highly unlikely in the typical usage scenario that your AI would generate more traffic. If you need it, you can [register for an API key](https://www.ncbi.nlm.nih.gov/account/) which will give you 10 requests / second. Read about [this on NCBI pages](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen).

## Usage with Claude Desktop | Claude Desktopでの使用方法

### Step 1: Locate the Claude Desktop configuration file | Claude Desktopの設定ファイルを見つける

Find and open the `claude_desktop_config.json` file:

- **Windows**: Usually located at `C:\Users\YOUR_USERNAME\AppData\Roaming\Claude Desktop\claude_desktop_config.json`
- **Mac**: Usually located at `/Users/YOUR_USERNAME/Library/Application Support/Claude Desktop/claude_desktop_config.json`

You can open this file with any text editor like Notepad (Windows) or TextEdit (Mac).

> `claude_desktop_config.json`ファイルを見つけて開きます：
>
> - **Windows**: 通常は `C:\Users\YOUR_USERNAME\AppData\Roaming\Claude Desktop\claude_desktop_config.json` にあります
> - **Mac**: 通常は `/Users/YOUR_USERNAME/Library/Application Support/Claude Desktop/claude_desktop_config.json` にあります
>
> このファイルはメモ帳（Windows）やTextEdit（Mac）などのテキストエディタで開くことができます。

### Step 2: Add or modify the configuration | 設定を追加または変更する

**For Mac OS**

Copy and paste this configuration, making sure to replace `/path/to/your/mcp-simple-pubmed` with the actual path where you extracted or cloned the repository:

```json
{
  "mcpServers": {
    "simple-pubmed": {
      "command": "python",
      "args": ["-m", "mcp_simple_pubmed"],
      "cwd": "/path/to/your/mcp-simple-pubmed",
      "env": {
        "PUBMED_EMAIL": "your-email@example.com",
        "PUBMED_API_KEY": "your-api-key" 
      }
    }
  }
}
```

**For Windows**

Copy and paste this configuration, making sure to replace `C:\Path\To\Your\mcp-simple-pubmed` with the actual path where you extracted or cloned the repository. Remember to use double backslashes (`\\`) in the paths:

```json
{
  "mcpServers": {
    "simple-pubmed": {
      "command": "C:\\Users\\YOUR_USERNAME\\AppData\\Local\\Programs\\Python\\Python311\\python.exe",
      "args": [
        "-m",
        "mcp_simple_pubmed"
      ],
      "cwd": "C:\\Path\\To\\Your\\mcp-simple-pubmed",
      "env": {
        "PUBMED_EMAIL": "your-email@example.com",
        "PUBMED_API_KEY": "your-api-key" 
      }
    }
  }
}
```

**Important Notes:**

- Replace `your-email@example.com` with your actual email address (required by NCBI)
- If you have an NCBI API key, replace `your-api-key`; otherwise, you can remove this line
- The `cwd` parameter specifies the directory where your repository is located
- For Windows, you may need to adjust the Python path if your installation is different

> **重要なポイント：**
>
> - `your-email@example.com` を実際のメールアドレスに置き換えます（NCBIに必要）
> - NCBI APIキーをお持ちの場合は、`your-api-key`を置き換えます。お持ちでない場合はこの行を削除できます
> - `cwd`パラメータは、リポジトリが配置されているディレクトリを指定します
> - Windowsの場合、Pythonのインストール場所が異なる場合はパスを調整する必要があります

### Step 3: Save the configuration and restart Claude Desktop | 設定を保存しClaude Desktopを再起動する

1. Save the configuration file
2. Close Claude Desktop if it's running
3. Restart Claude Desktop

> 1. 設定ファイルを保存します
> 2. Claude Desktopが実行中の場合は閉じます
> 3. Claude Desktopを再起動します

## Acknowledgements | 謝辞

This project is a fork of [andybrandt/mcp-simple-pubmed](https://github.com/andybrandt/mcp-simple-pubmed). We would like to express our gratitude to the original author for creating and sharing this useful tool. The enhancements built upon this solid foundation aim to further improve the functionality while respecting the original design philosophy.

> このプロジェクトは[andybrandt/mcp-simple-pubmed](https://github.com/andybrandt/mcp-simple-pubmed)のフォークです。この有用なツールを作成・共有してくれた元の作者に感謝いたします。この堅固な基盤の上に構築された機能強化は、元のデザイン哲学を尊重しながら、機能をさらに向上させることを目的としています。

## Verification of Improvements | 改良の確認方法

After installation, you can verify the enhanced functionality by using these search and retrieval features:

1. **Date range searching**
   ```
   search_pubmed(query="cancer therapy 2022:2024[Date - Publication]", max_results=5)
   ```

2. **Full text retrieval with text formatting option**
   ```
   get_paper_fulltext(pmid="39661433", format_as_text=true)
   ```

3. **Enhanced metadata search**
   ```
   search_pubmed(query="BRCA1[MeSH Terms] AND breast cancer", max_results=3)
   ```

> インストール後、以下の検索と取得機能を使用して改良された機能を確認できます：
>
> 1. **日付範囲検索**
>    ```
>    search_pubmed(query="cancer therapy 2022:2024[Date - Publication]", max_results=5)
>    ```
>
> 2. **全文取得とテキスト形式オプション**
>    ```
>    get_paper_fulltext(pmid="39661433", format_as_text=true)
>    ```
>
> 3. **強化されたメタデータ検索**
>    ```
>    search_pubmed(query="BRCA1[MeSH Terms] AND breast cancer", max_results=3)
>    ```

## Last Updated | 最終更新日

March 23, 2025 | 2025年3月23日

## License | ライセンス

MIT License

This project is licensed under the MIT License - see the original project for details.

> このプロジェクトはMITライセンスの下で提供されています - 詳細は元のプロジェクトをご覧ください。
