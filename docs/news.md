---
layout: page
title: News
permalink: /news/
nav_order: 3

---

# News

{% for post in site.posts %}
* [{{ post.title }}]({{ post.url }})
{% endfor %}

