---
title: "Covid-19 and Output in Japan: Chiba"
keywords: sample homepage
tags: [chiba]
sidebar: home_sidebar
permalink: chiba.html
summary:
---

{% assign fig_loc = "./archives/20210126/Figures/Chiba/" %}

## Last update on January 26, 2021

Replications files are available [here](https://github.com/Covid19OutputJapan/Covid19OutputJapan.github.io/tree/main/_archives/).

Link to other Chiba pages:
<table>
<tr>
{% assign cnt = 0 %}
{% for page1 in site.pages %}
    {% for tag1 in page1.tags %}
        {% if tag1 == "chiba" and page1.name != page.name %}
            <td><a href="{{page1.url | remove: "/" }}">{{page1.permalink}}</a></td>
            {% assign cnt = cnt | plus:1 %}
        {% endif %}
<!--
        {% if cnt == 1 %}
            <td>here</td>
            {% assign cnt = 0 %}
        {% endif %}
-->
    {% endfor %}
{% endfor %}
</tr>
</table>