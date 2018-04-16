# Metagenomics

O processo de construção e de elaboração para o pipeline para análises de metagenoma se baseou no modelo do [EBI](https://www.ebi.ac.uk/metagenomics/pipelines/2.0) (European Bioinformatics Institute) Metagenomics 2.0.	Para isso o pipeline foi dividido em três etapas, controle de qualidade, análise funcinal e análise taxonomica.
O fluxograma do pipeline pode ser encontrado [aqui](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=https%3A%2F%2Fwww.draw.io%2FMetagenomics&layers=1&nav=1&title=Metagenomics#R7V1Zk6O4sv41FTHnoRxmh8fumuk5E9Ez03d6Ttx7ngiMZZtuFjfgWubXXwkQBimNsRHIdrnO0jZmEcr8clNm6kF7il5%2FTb3t5vdkicIHdb58fdB%2BflBVRZmb%2BB9y5K08ojq6Wh5Zp8GyOmt%2F4GvwD6oOzquju2CJstaJeZKEebBtH%2FSTOEZ%2B3jrmpWny0j5tlYTtp269NeIOfPW9kD%2F6v8Ey35RHbdXaH%2F83CtYb%2BmTFdMpfFp7%2FfZ0mu7h63oOqrYq%2F8ufIo%2FeqXjTbeMvkpXFI%2B%2BVBe0qTJC8%2FRa9PKCSTS6etvO7TgV%2FrcacozvtcoNKBZ%2FkbfXm0xHNRfU3SfJOsk9gLf9kf%2FVi8ICK3mONv33bRtiKhQr5v8igkn%2FFH9Brk%2F0dOmxnVt%2F9Wv3xDef5WXeXt8gQf2j%2Frc5Jsq%2FOyPE2%2B10QgR%2FhXrN46S3apX72DXh4ir9I4p5qDX1ESoTx9wyekKPTy4LlNdq%2FinnV9Xn3plyTAT1XnFatrdkXHitE1Sld6i9xL1yivrtrTAX9oDGN%2FqKAOTKnqjZ69cFeNNkM%2FstnKy%2FIfHAn3BCIz9rIJcvR16xWT84Ih2ybSKgjDpyRM0uJabekhe%2BXXU9%2F4xfRttFjhX7zUr0hnzLso8ozSHL12EuCVmbhqIhWz%2Bv6yh59Cz9k0oGfOD9OsNdsdU6sagzHQ5PkYP75kel116IH%2FFgfmyim8X1z4BaUBfiOUVk%2Fqxf7VC5XcVx6yRkHEhzT13honbAmnZ4cBoyhtOqvWnKFUecdzUWJwKMEfg6WXIzfyttsgXruY3dFs%2B8aRHHNq3qZkirLgH29RnECmvno5fLbx8cH4mQABUy%2BrRB%2F5GgbrGH8O0YrcirB%2FgJXKh%2BpwTij7McNAxOP4uyDzow4A0ED2UocAaKsLzTQFIc5oE0IBAGcAgGMl3DmAszgqJXjgYbCYfcO8Y4Zk8pbBM%2F64Jh8zL9qG6BE%2FxsOU9GZ59uz6SZpizY%2BWsxyTrbpmkdJL6BE8ksaN%2Bt17VrCA2FuGybr%2FHQUKc2x5qD4ozJfmwjTMljBXDDGsxchyDWCtmt2arOUIYC3F0EcS5kZLlI9vwti8DNedixDijtGWHYouVohjTT%2BqUXqApo9YP5vWIb07PbUN0cTuiyGbNzVZkfaaX7D%2BPG7AChByBmOwaj0NViH6U%2BFmX4iMaxurxglc34vHFerKN5hc0WRxOR1Nk823YZC72A5JvTRAmVs4V3dbUTUZq12xAGbXxmJ2sQq9jkiobEii6bg5KqMZFNMeqhn6YUQDMCJc7ff2jK3hrnGP8JDSoAQlTIfhJUM%2FQ3SpWXN6wuiONg4qmqCYd9ChnEovzT%2BQcCs%2BECcxosc%2BBWTg1fVLeoYfelkW%2BOXB6pRTwiGrJM4bAm4%2Btz59%2BjSEfgTlspSPximfTZDlyTr1oqzp0jb9yiLeF3sDXFN6COuVuL5tQ%2Bm9udhRbT4e37F57qX6sAJ0nGWrbR1nTBmBVIYruTM8HmnC1AGEqURL0OHA6OOpzd0CcXfzz3SUFjRUwPxzxrL%2BhocCesXmlRnR5lNYeKoJWHiKNAvP5Hg%2F2pGF1hC53zATu1svSNHSxSr7wqDgJxExJopn1Muo5Mua2BnVZyl4qcFQ4UWHVIk1EmBssbGBgwu4IiKjvfCiA6Exy5ZmePPBsW3gf3eTLYrdFK0QfhcfuUm%2BA%2BEyaGH2OPPuY%2FmC1mUVhYn1QrJfsUcyi%2ByRXH95vAwYPnTyZDiRMC%2F7YZJhkX8KN8%2Bvgpttk2Fm%2Bn0SZp7AY5%2BUmalWawlmaaEqOhqWmZfIjZPn5HYkst5HIutjOar2WLky4%2FKqBRgRii2NWfX5KJqtHUe12DiqZrNzDHo71xTWg%2BhqS1s41Q2x3mot463DYdkWSQsUXTv9dEOaN2zx1v0qwK9LtAdJIHMjXz0QnQ0wDV%2BHpg7F6IWYXQIivWTAOfHCyZDdF%2Fz5dbYIkkhA8JgQ%2FIHkXlOamz92SXkCoX6R7Lw%2FVF67fYu9LHcLVx%2BbluQFG88t73hqiJkbYIq2%2BM75LE%2BRpNfMcrTVC1OjcwDQ%2Bx4dV4sRB2a9hcnaneO%2FgbdpMdlQ9gI4FltubsU5Ky8IdynKhj6F8oiYdZRzeURp8Yg4BAwZ184PdwSiWBkUGMXzn8RJ9DZ4lNMuBY0QAqHRWZrPCHiNdR1FK5yniEhotK3BNkWvCg0uH%2BL8RaCG0Wlxdn0%2FuwAIlIxgr5%2BXw65TV75iCNNkimzYC1TH7LqgnTDJX44h27rcMhi%2BElhlYg7ntmEWrHV1JizAqjJNWD6oB2l7SPHevIE4ZICV4XCLVusRxS%2FOhBVgnpCYHh6jEPPp2g0TXW%2FrBQ2KBELrjIqIhUbFHB7Q7mWYWBPXjopfkDmg1w2N0etaP71%2Bcr2G044Ym2xR8LGBnXwBfRNRFSGKbZ%2FBa2cVJjdsV1VgyU%2Bb70%2Fiw0u0gdkSoKMmsDLMBMYsNlPM%2FZ%2FFGsQzrfFnjGUe2%2FwiU7GuhHU3AvSd6DCLTAOt0r8pWt118J4xbauNA1AFm4AKFlHrKDg%2FQuiCmw2UoEtdcDtHgZy1bnlNvitEJX0uvFVA7xQJvny%2FlHjzbBdFXopn47EQhBwt31lGrqLbB0zHKSqy8LsMBtPJVUHMYvbcaQNv%2F%2FvEuewOVBikSctlVw1zVNpMWw%2FfzpEWnvQJG5ua3o73PzosaEr6izAoHb4OqCg69YdZRm07R04SNCOhHh0bEFFjFdTUzNKY2C9FJvsjVrP85G6SaIGt2ZNNTttHsMm5sA1dVHOkOqebTqXddyptEdLemGgB6pRERilVS5BUqvvoSChbqjlpz%2BN%2F4xnBR371MCPyVtL1yRCG8S19QhmiGuMGOCl7Ky3mnqmKrB4kigoVn9MzJUecDGYR9NGq2kgJ6zojoHPQGeQ%2BpcvGTdLVYuhq692hRO4CSxEcbFZ5o%2Bwr%2BvElRdsrMxzYBaKJDQdzeCnaeZkrl2c5QAXPmjkxAAf4KfUbNDDxqWw9Ov%2FqpwjFV29ssCGVR7rqNY2xYd58B4AW85vCIySHnHmGqrVMFO%2FMKxpvkP%2BFYi%2B6QUu8Nscm8eZN3pt%2FSiIvXibZ1U%2BtwjaRflS0KSVP3bu8FSrJcq9rak%2Bv95ymiwswmeqUraTrdJvGZH5If%2ByC505OPaN8dpImh9B0Tgp82v6jMZ28moQW0xozB%2BQnH84taihKc98QhLYHVesjZ7XlhlJFe1erNeYb6stMjw1OSjLbFDdYV%2BCAvuTvxHYNNHqmN52jeS3eZVsiP8genrSHjwYSLMomcd4A7ClThr80VWh2wbgF68N6w9kALjUqaSQEdm1eI9P8mS0J61Zpu0mK3CV%2BcJoF%2BZvr4Rl5y1BZ5F4n2zTOB27BJATtKU2zf8gPj%2BUi%2BQd8grV95TOD2jcpj%2F1EgEdeY%2B7h%2F33bxeT7B6v4fy3Bh5Z4DPMD2UN5Wpz3AaN27ldm3r%2BO5x633kukqp3GDFT0doacBbVNn48EdkMdHqs5rZ8kG%2FW%2BtrIPhWZwtYSGJk1oGJBlJLxP7nsgobT0p3o83XI%2F8r4jF0VblCbpKcL%2Bp5dCaKDlv%2FpfsouPXjRE0l6A%2B63r7QxOx%2BlpY4kIkRs0rDodaK%2B9jy%2BUb2XI2ytBATJ4KFKanXO9cLvx3BS%2F3Mrz8yCJ28iFmuwe6dt7qu1mQ7bbT38Wk5LtwtxblnZZlqG93UWMNGKLHbPkkh3G6In2GvQyvdPDL9WMs5l9ISa14jRrFHHSlYjZNgFOWcweJAbqKWypbvGLCr1n3u6e%2BUqqHpfc3FS3hHOCpQwjnJW2cK7OkCCZa5KMXjfH7b0yt5yZMzcdy7R0XbHN9h3FhZnqVwSEfUNsFQbaMsBUIs0TF8nrNkxyxi%2B%2FbimnqBZLAqDH7WhyTp3zO8jx8NuDokJfr5jwacLufJulfzHMRIHfOYMqe35u4FeZrl8FtJHCH7sIu0bUMkmRt%2BTDv236qwxzbLwtOS96XZP9jWerMHnxN1jkzjCotyFJLQJRyYaJDfIfCIplDaMYKJrMvgWaDvQ0pdVP4rMWxt0fRt7ue6oClEapqry6DnPcuo6D9sfF0UB8GUj%2FXTr4%2BrTf%2Fvjt6bc%2FO1T4SUkIwtaOGFn%2BuN%2BhfJpkpnFLBgR4H2ezJNCw0pEnFQQ0le4pf%2FGTLiCTEpx%2FiQLhvK0yqYGxTHxsq5DXPSYiFiXJPi%2FCY7N1guXAyoj5lBVaqjaRjJi8UFG8OIDtbC61UZmPl9uoAknxPQOVJSKy2arMFxa7sVgNJYy1LCdFT5eZSaU6dg%2B00WQg4an4Gk%2B9%2Bexnb5m4H1OCAFH%2BEaZCT%2BdoknRArgACLEKF%2BnOIKSCbpE5%2BPnPmreWeetGdVc4yDScNagii9o2BiFfdGm%2FLKzNSPumtUXSreKhLcjprHEbDgwlE6cbAQ9Uc4igeRud5aKs3WiQgg%2Bf5zd7U2VMS514UxF6cA3mSN8L4OtOgDq5fGI3xteHrcv0YX2kzvmMf1QT0mr02cOptQibHy6X1Uhm3gHEUF6WReqfrDDMc4oZ%2B1NFVSINfRhUx17PFsrqriI9dMLycnE5Xc2WCYyd5fasECFWNFao2YE1A8QMhTaj0cY1riaayDpbryxODOlAPf2OczHXYAZezx%2BJkPt3%2FxiSF7siUFMDWgDfGvwbX6GFK%2FjV5%2BeCHyIvdMvJHA3%2FiPIuFlx2K88kJ7rGVWnrfIkkhwT2aiyOvHSOvI6X4Dya0y6M8d9uaqP9JV0Fde5d7act2EGVU4f1P%2BkusG%2BoOWDcLeKu%2FzxTL1K25YduqZtCgzzSiiHd6bkzValzHb35%2BR9O0QJnzrVmKLDdPachYsjpddkfoJs%2B7snReWuvysiysaZbRePummyz81nZS6CJvCU1Atcl5dmc3XURu23I%2BXajckkEXndMSZMuOasMOQW4Y6a4RJHE%2FN2yijhkKs9GRAtVyjlaUAGwOw%2BPj9F41HXqg56ZCZ7HlNFUHDlvK0zODCehbYx%2B5k6DNtJS62JR5kLCFA5rQPxYftZ31QlqewUi97YGpyleO0b9%2F%2BQqzrFHHV0RzEjtkmv4lipMMZXh1Rr8d%2BbqaG19TdTm4m4fuyDOxOqrLj2RdLlDuNboD5RtM1%2FXGBUoR%2B9Rhj1KRDnYTElKRjq8gxejLJHjHhejtTaotG9i1fK7xApmVdmf1tRCwQ9Cp1Um8O3ZdkgfyIeT1oqnHc0a%2BN92ri%2Bwk7A2VPNcHPCZA6NBEkEk6QMz5FL%2By9QieVNdPQswmePY4J%2FAliEKvAAjhYgpTMsX%2BJgiXn723ZEeGieHjf6ffPm4STGaSPrg3OzG6KqioJuAH0mS%2F%2FZVfyR2pP0PilOgLnX2FOfS799o68bOX5XSU5M22WVDGOMmFmAXXQfwxyfMkovKgestP3KDIjpeFt7ssJAE5eU2kQX33iMiFSqiMmKbItnzUadVA04HVAd5R5roI5uHL6rG3nyhgEVMVZm7PRyVlm2SvDp0SUv5cnPOzvj%2FyV%2FWq5FCCryWBCPx1EyyXKC5UR%2B7lByPcH%2FGMPRGFYeCBP%2BHvyv57Ef7eYoH%2BlMT4XbygIB7CrPWCshwk6x5kxwlLTXtgEQaio2qKICO%2FupUsMpQ%2BV5u23sk5lJwGL9JHJCefHvTlZ%2Fdlk4RYw6YI3ck5lJxWTykrhJwKHybkmoOVxtJdTV%2BqmmbizKYGyAMTqOWjjNYR1%2FoL%2BbkXr0PUeBxTqWtZ0ON68KsX5iiNsQX4kVjRQ%2FeE0WhqTYOTOR6eFeb4XTrV0qmG%2FyDbAWIuMdKJDz57eKjeGh2QSe%2BFlP29%2BzNIDNkTo5GYri9M19f9gQsKTxJTqVu4t5oEGtISGQxleH7JiRugXXuLZpCE9bKuhD6PKi8fmbC3X7%2FrPrRN3poYL1y02yMK04281yNRaiDkds7jFiQKhJ7Jbm3HguL9wv%2FdgyB%2FwCCCJHILKe2WocK30oORN6IAA%2By1slTkjSJM1i7%2BYarJWK3gYZSUIUwyI5Sarf%2BRP5rIn3Ish9BThLMrK%2BikUQyJYffydg7JWDl7cClGeycLqNhHsaEopph9a8bpqf7QXHuyWotPkJadSD0CHQZq%2B1OCehRQinJ07tVLnnsafJJhmvBLP4NME2DZaAob5WAwrKNDlCxdfY6yfl%2FKwKQLU9TTBcolFFqi0gqVaQIcLnUCcXRTaQh1Sf%2Bl%2BFsaHy2%2FHCif%2F9Ak37mMoS1ausIP3m%2Fo4%2B7iYJV6vruMmq9%2BGfIVGObWl%2B0mHpm5A0M63tJw3EGeNW%2FvS0kZrMMC7LNZt5ltOSxskvF5Mo6PKb0wRKy2%2B8LETEoXlCPRfW2w99rgId4bt%2FUesxuuRbuqNHlMA5zirjXD%2FizG5xFTlkrRjx0pfcYfSkvgvsbTY42nxuyQZTyQ3kLWeDQ%2Be7MZ5Hyn9BVCR2CtbkQ68hl5gGl01w63ph1sU%2Be5TIey9YVoBz5h8K4dBmkHqzcDdGgHiN5ipAqfUXjXDmLoCGmH0ehojFPf93AwvCUv2K5XlGinckjr7VR3fZkuieb8NiAjzLwlbcPdejwnL3MMDv8Jj%2F19VvtE%2FkYNF7UGIzm%2BNpvNLmQqzEuiiymCLkXuQpVuw%2BUwvOPQn22316cURQUsb6hcr95UcZgK53OoIVLd3bmrcudMqqRqrqImd3PZ0zmxQqA%2FU%2FEBZfw2nlu07rmn5LNWvtGfvh0p%2BRA1xVj5fOz2Tk2R1ISy70ejpj1KSsJpzXmncRwMyHEw5TkOdndPigoKZ7jHrQSP5Jlo11aCR%2BnHLZkzZGR3GML7SVSXsq2NVLavV73%2FrPg9wOu3Or1VBTYUyJS5da5b2S3neJeKk%2FcpvPy%2BFUrd%2BaXLaBmvcwXdguEdCEagGbZmSmv1Uo%2FndPys8fxv3QxbCcEKGwyxj%2B7QKQ8AFcHjQYfepGfHwErTNWftAD%2Bfh6Ke3QSPNl5T7L6N16bqMciQ2dYY6vXvMci0ubbY7uvCegwyO6%2FQWoyDI2MvsAS3t9RMPuIxJrOqo3Ar0Kzr8riVaSbkKGdzK9MR02al1ljc6lSmnURuPdyAsaP14Umh3k8BCpfuzwF2GbACzWbb5ep4tBd6%2BGntGc8Y41fsx2f9otGDmkqeM8KftvGaxA8wUp%2BrD3gie3R%2FfMfBb4Xdl0uZA9stKJCQEmO38HswD2fSRfJKguY9YCSI88rHZv2xcee4I07miBwndvG4JH3hBXWS%2FlKIPM12aIrBWLegOwRuQiSEyI6gQN9JyRhXVeVF0zdaEQjqgsnIphGrCD7%2FsvqKLg%2BMAqDl0OJ2iqy6YdM0yOourT8fWVeFHihWoVGXUAJ6OprEn4OeVZJGXu4G8XaXt8N5F4IjKTE%2BDnn6fELk1RuKSmusYF1d6ygoSFN3bpwepvp8lAWO%2Ftm6N0JAeQ02OlrinyNn013shmiVobuQPShkjSmFbC0cBBGYhCFI3cydvpS%2BlilTidayQyR9sYhcJuvUi%2B5UPkjlSVGsK6Po2dOWQK9Jz%2BpzIFFLl9dETAe2oRgC0y%2B%2FPf31n68XGJwTADVFt9utyO0pI2060Mt%2BCKW%2B%2Fv3h9y%2B3SSeTSYiblk6O2ELBWxF8qgoIPnmBnHo8guCkG3pvMG2SaLHLprVGptkV2GB2T3jUaCzl2K7AbE7OeTQVa3Pif72I0CReZOSfP%2FpSuKoIaJCyKMo5mJZPMu8J8rJqd%2B6H06oE%2Fi4g%2BqiPK1bNOUNbFejcBDbnEEJasXbK1%2FdNSp0lpcWTElKQQkhJF6HH9hmuSj1qCqAeNWnbJ9TjAfC2TZGIVI4lipMgQy7Re1uUdrn09eHWo9%2BZo6%2B3nY%2FHuqt3q7srUEolxKrV747%2BiYAGCg90XZ6jD%2FQdEgzoIF6R1tVuBezUTXb5kWXNO7BV1XGU48BW6caE4oGtjwHs283nqTHcBvZU%2BTz%2FwcD6c%2FEN%2BQSeobdA4ShY3gb%2B93L%2F0rPQGwbx93JcmzzfZsUW9J%2Fwf9dBvtktZj7p9%2FEpJlvMR9jsxe%2Bgfvod5d4aSw6M0gx%2FXYTJgry6l%2BUoxR9S5C0jNIuW5QP27DtvzHPF21cbuFKZKgJYxwuKXOGvDWY6JB1Gab10w2pfB%2FKVdF1eIoR%2BOF9JrB2fuvh%2B2zTxUXam0Hg%2FKp8r2X6ESgGgXXvEqPzDS%2B%2BCeeLOCCdmb4OMMNpShd6dvC3E9qsF%2Bk1Id2j1VlXlSffDUVFBSM7TXewTr26FLTHP%2FbHzQhdDpTvR5g7tRqPf2paDYqzKWPGakdqVXnUyuW5AQVZD3hok0JBuCFb%2FQlHyjOc1wRf%2Bzy6IUOod3j1RGhpFqE2TcZS0uTGh2jRG2lr4ZkMmFGJt3ElLLq7H009t%2BslSTEzFS70wRKEbYPM4D1Zvrr8hIA18N0M%2F%2BilUaCzvW88qJpNzoAMpBzXuhcsCepN7gKSvLAB28tSNvu08RpAFYnNGnlqQvoY65YmK0S0mfcRQgFp0MOQhgsqmOr7Kbnu6jQqv67SVHQCnlBYycOrI0NkYKzlKS%2Bf3rqDPQ77DKOg65jWFsW5QvhGNfKvDWlev21w3ITfZlNaxQAe7mI0O%2FRRtsRLP3SJ9Dy0rY91bpQjNCnlwW8JgIjOAFQagtT5WSxrdvC9nnigKoOVMeX2l9Y4WceNbAYUkiPB83S2Bs8C%2Fz%2FzvsgTGAz%2Bw7Kl9uhwiTCOBVdWRaI7pzmFdXrRqa1CCopCk4D9WSfQfCHTxy5OcIKBcg%2BvDuEC5dxU9k2mRQPVyoy9OMU3DjBETjQAmUEdngkLM9rhPWrzawdtU2mLl%2Be3b%2FBuFz4iUcfCXestvuywn8R83SQM8TKwCk%2Fj05dJb5T2LGmF0fxdIBwC7OIrhvAvYSvCqDECoL0xdHCLBAOzoC9MU5l60DQM%2FiYs0pC6ZfwN4omCpi7mcSeMr2hiIumGXyoFyiBx5lV7O4RyiFqIyYka5ufeaxEnUuW3L9WOq3hqO2kcOVMc8HqZGqcm4ZUwBxVYGdSRlYIovtioYuLHQ36osvwXM6CazD%2BO0mFFH0UPzDsxUSwBXu8JHG8y2QDOXCJrDK3ycS3%2FC9mGHO8%2Fz3iH54bEs3SfOob195f27n5YI%2B2kkJ3Y55Y4HlwFym2k3ZkNxA6i%2Fhy4C48o4%2Beocxq8JxjViW507NWkwrsfTgPE2DWI%2F2Hqh6ydJugxiDJ5SDV4%2FIqy5MVNbmLBsYMdazeIxQYv1h%2Bm9fg5xtYSREccYZe4qTSK33vf8to35uuM9NUwsoHHNaBvGGRSKE6YeXbVdYqjAhgmGKm3Lxno8DYRFuzAPMJDcFL%2FLyvNJhPdWBJrNLNKZqsqLMx1AiyZgy2aDhhbvaOmLFiCcZKjSwkkG0KgL2vH3BnCitHFi0K1dj7m7QtS%2BNtJChtoXKNe2klFjogUU2hlUBlD4GJGfhEUHnQIwNwqTegnuKEyEGF%2F3MkgeCZTpWx6jKq0cqx5P08DyvreMK5f0f78VC8tU2it2hqIDFhZQGUxjnINmWx%2FFH7ndxQVDgyws2sRSBlx4CyvbRZGX4okg63O3ojksGtqoN0zVplxQGKfj4e2mitSQaKkVXZ7fTsfTDESSXUQIRNwSMbfiizi0mK4TKlBzZjFQ6ZWVXej0RYL1eF60fF26OamlaMUehWVHXwZdVNrhod5C0uaDKaP1djT0iazfmxFhOpDtZkzW2xEgYK9GBd88%2F%2Fv3OFhhRHUtjt4atnQmRwdsYj5Wgx3DuJsHJ2ILKPw3JJoHOl%2F43%2B56ev0QYettFMsGPM3x0tjUMSByw54m2F7OkQgR3tMsIFKVxt4qSiZO9hylm1QbJWoXTK5ut%2BE6s7Od7imtKrVe%2BRmThFYXCdUbIaG8GLTBZy2BVb%2FXL%2B7Y3fGmFndANgUJDGw3b2GyRvGthGRslUmHtoGsoxGnmbdui2lu53zdwjQzkS9nxGnGX9OEJAPXv%2F2Kp2Xze7Ik2zD88v8%3D). 

## Instalação
O processo de instalação do ambiente pode ser feito utilizando o Dockerfile ou pelo arquivo de instalação. Para mais informações de quais programas são instalados e suas dependencias veja [programas](https://github.com/nandomaciel/Metagenomics/blob/master/doc/Programas.md).
 * Shell Script:
    * [sh MetaInstall/Meta_Install.sh](https://github.com/nandomaciel/Metagenomics/blob/master/MetaInstall/Meta_Install.sh)
 * Docker:
    * [docker build -t Metagenomics:0.1 MetaInstall/ .](https://github.com/nandomaciel/Metagenomics/blob/master/MetaInstall/Dockerfile)

	**O processo de instalação e de configuração do pipeline, dependendo da taxa de conexão pode ser extremamente longo.**

## Controle de Qualidade _(QC)_
Para o controle de qualidade é necessário verificar e validar as amostras de entrada. Removendo os adaptadores, sequências de baixa  qualidade e a remoção de contaminantes.
Nessa primeira etapa do pipeline é utilizado alguns programas que facilitam esse processo entre eles:
* [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): Para análise da qualidade inicial das amostras informando quais sequências têm baixa qualidade entre outras análises.
* [SeqPrep](https://github.com/jstjohn/SeqPrep): Caso as sequências sejam paired-end e caso seja necessário, é utilizado o SeqPrep para mesclar as sequências que se sobrepõem a uma única leitura também pode ser utilizado para a remoção dos adaptadores caso solicitado.
* [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/): Para  remover os adaptadores incluídos na sequência.
* [Fastq Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) : Para remoção de contaminantes indesejados na sequência

Após a finalização desse processo o usuário tem como saída informações sobre a qualidade inicial das sequência, dos contaminantes encontrados em cada amostra e informações de qualidade depois de concluído todos os processos suas amostras com o formato:

 <p align='center'>SAMPLE_S2_L001_R1_001.fastq.gz</p>

Onde:  

* SAMPLE: Código utilizado pelo pipeline para identificar a amostra  
* S2: Sequência do código de barras ou identificador do código  
* L[0-9]: Lane do sequenciador  
* R[1-2]: Número de reads da amostra paired-end ou single-end  
* 001: Constante  
* .fastq.gz: formato do arquivo zipado  

## Análise Taxonômica (16S/ WGS) _(AT)_
Para realizar a análise taxonômina para amostras de 16S pode ser utilizado dois programas o [QIIME 1](https://github.com/nandomaciel/Metagenomics/blob/master/doc/QIIME1.md) ou [QIIME 2](https://github.com/nandomaciel/Metagenomics/blob/master/doc/QIIME2.md).
Para análises 16S ou 454 pode ser utilizado tanto o QIIME1 ou QIIME2 para outras plataformas como Ion Torrent recomendamos a utilização do QIIME 1.
 
## Análise Funcional _(AF)_
Na análise funcional foi utilizado o [Kaiju](https://github.com/nandomaciel/Metagenomics/blob/master/doc/kaiju.md), [FragGeneScan](https://github.com/nandomaciel/Metagenomics/blob/master/doc/FragGeneScan.md) e o [InterProtScan](https://github.com/nandomaciel/Metagenomics/blob/master/doc/InterProtScan.md).
 

## Sample-Metadata
 O sample-metadata é um arquivo informando quais são as amostras que vão ser utilizadas no processo de execução durante todo pipeline, para criar o arquivo [sample-metadata](https://github.com/nandomaciel/Metagenomics/blob/master/doc/sample-metadata.md "Como criar o arquivo sample-metadata").
 
## Amostras 
Para execução de testes do pipeline foi utilizado algumas amostras 16S e WGS retiradas do [NCBI](https://www.ncbi.nlm.nih.gov/), tanto paired-end quanto single-end. 
 
SRA/16S | Description | Layout |
-----------|------------|-------------|
[SRR5028922](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028922) | Wheelchair in Shang De Long-Term Care Facility | PAIRED |
[SRR5028931](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028931) | Computer in National Taiwan University Hospital Hsinchu Branch | PAIRED |
[SRR5028971](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028971) | Bed rail in Nantou Christian Long-Term Care Facility |PAIRED | 
[SRR5029132](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5029132) | Stethoscope in National Taiwan University Hospital Hsinchu Branch | PAIRED |
[SRR5028967](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028967) | Waed deessing carriage in Nantou Christian Long-Term Care Facility | PAIRED |
[SRR5028961](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028961) | Anamnesis in National Taiwan University Hospital Hsinchu Branch | PAIRED |
[SRR5028953](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028953) | Doorknob in National Taiwan University Hospital Hsinchu Branch | PAIRED |
[SRR5028945](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028945) | Sphygmomanometer in Nantou Christian Hospital | PAIRED |
[SRR5028934](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028934) | Floor in National Taiwan University Hospital Hsinchu Branch | PAIRED |
[SRR5028965](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5028965) | Dialysis machine in National Taiwan University Hospital Hsinchu Branch | PAIRED |

SRA/WGS | Description | Layout |
-----------|--------------|----------|
[ERR947518](https://www.ncbi.nlm.nih.gov/sra/?term=ERR947518) | Metagenome analysis of Amlakhadi canal | SINGLE |
[ERR947519](https://www.ncbi.nlm.nih.gov/sra/?term=ERR947519) | Metagenome analysis of Amlakhadi canal | SINGLE |
[ERR947520](https://www.ncbi.nlm.nih.gov/sra/?term=ERR947520) | Metagenome analysis of Amlakhadi canal | SINGLE |
[ERR947521](https://www.ncbi.nlm.nih.gov/sra/?term=ERR947521) | Metagenome analysis of Amlakhadi canal | SINGLE |
[SRR6037376](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6037376) | Diagnosis of respiratory disease complex in poultry | SINGLE |
[SRR3738816](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3738816) | Discovery Deep brine seawater interface in the Red Sea | SINGLE |
[SRR3738821](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3738821) | Discovery Deep brine seawater interface in the Red Sea | SINGLE |
[SRR2657586](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2657586) | Arabian Sea oxygen minimum zone Raw sequence reads | SINGLE |
[SRR2657589](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2657589) | Arabian Sea oxygen minimum zone Raw sequence reads | SINGLE |
