class Tagging
  include DataMapper::Resource
  belongs_to :item
  belongs_to :topic
  
  property :item_id,  Integer, :key =>true
  property :topic_id, Integer, :key =>true
end