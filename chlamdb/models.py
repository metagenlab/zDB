from django.db import models
from smart_selects.db_fields import GroupedForeignKey
from smart_selects.db_fields import ChainedForeignKey




class Database(models.Model):
    db_name = models.CharField(max_length=100)
    def __str__(self):
        return self.db_name


class Genome(models.Model):
    genome_name = models.CharField(max_length=100)
    database = models.ForeignKey(Database)
    def __str__(self):
        return self.genome_name



"""
class GenDB(models.Model):
    database = models.ForeignKey(Database)
    genome = ChainedForeignKey(Genome, chained_field="database", chained_model_field="genome", auto_choose=True)
    genome_name = models.CharField(max_length=100)
    database_name = models.CharField(max_length=100)

    def __unicode__(self):
        return '%s %s' % (self.database, self.genome)
"""

class GenDB(models.Model):
    database = models.ForeignKey(Database)
    ref_genome = GroupedForeignKey(Genome, "database", related_name="ref_genome")
    query_genome = GroupedForeignKey(Genome, "database", related_name="query_genome")
    genome_name = models.CharField(max_length=100)
    database_name = models.CharField(max_length=100)

    def __unicode__(self):
        return '%s %s' % (self.database, self.genome)
