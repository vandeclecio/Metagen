# Sample-Metadata

O arquivo sample-meta é um arquivo .txt ou .tsv cujo objetivo é auxiliar a identificação das amostras no pipeline na etapa de controle de qualidade.
Na análise funcional o arquivo sample-meta tem uma grande importância pois as informações que estão no arquivo serão utilizadas pelo QIIME para gerar as informações necessárias para fazer as análises alfa, beta e das OTUs.
Na análise taxonomica o arquivo sample-meta é utilizado somente para a identificação das amostras no pipeline.

## Construção

Para criar o arquivo sample-meta a primeira identificar os IDs das amostras a segunda coluna é informado os barcodes e na terceira os Primers utilizados nas amostras apos essas três colunas iniciais pode ser colocada qualquer informação dia, mes e ano da coleta das amostras, tipo da amostra e a ultima coluna deve ser a descrição da amostras. A primeira linha do arquivo deve ficar como a tabela abaixo

|#SampleID|BarcodeSequence|LinkerPrimerSequence|SampleName|Treatment|Year|Month|Day|Description|
|---------|---------------|--------------------|----------|---------|----|-----|---|-----------|

Apos a primeira linha qualquer linha que começar com **#** é considerado como comentario.  
Exemplo de um arquivo [sample-meta.tsv](sample-table.tsv) pronto.
