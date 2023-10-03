use crate::field::F12;
use crate::point::PointExt;
use crate::parma::{RBITS,RBITSZERO};
pub fn genpairing() -> (PointExt,PointExt){
    let xa = F12::from_i128(49015138457);
    let ya = F12::from_i128(36336695411);
    let a = PointExt{x : xa, y : ya, z : F12::from_i128(1)};
    let xb = F12::from_i128_vec(vec![52649681422, 58330021724, 18361934063, 37544484019, 48378959847, 34344081429, 20873394922, 26784166228, 25856414820, 31733507710, 1579608448,25687127276]);
    let yb = F12::from_i128_vec(vec![36962872648, 28863545577, 49052618431, 14807993877, 27927561526, 25370717733, 16963613074, 45537344140, 11799967561, 30748944721, 10754770505,42033372585]);
    let b = PointExt{x : xb, y : yb, z : F12::from_i128(1)};
    (a,b)
}
pub fn miller_loop_t(p: PointExt, q: PointExt) -> F12 {
    let mut f: F12 = F12::from_i128(1);
    let mut g: F12 = F12::from_i128(1);
    let mut mp = p;
    let mut mq = q;

    let bits: Vec<i32> = vec![1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1];

    for i in 1..RBITS {
        let lmp2 = -(mp.x * mp.x * F12::from_i128(3) / (mp.y * F12::from_i128(2)) * (q.x - mp.x) + mp.y - q.y);
        
        mp = mp + mp;
        let l2mp =  q.x - mp.x;
        
        f = f * f * lmp2 / l2mp;
        if bits[i as usize] == 1 {
            if mp.x == p.x {
                let lmpp = q.x - p.x;                
                f = f * lmpp ;

            } else {
                let k = (mp.y - p.y) / (mp.x - p.x);
                let lmpp = k * (q.x - p.x) + p.y - q.y;
                mp = p + mp;
                let lmp = q.x - mp.x;
                f = f * lmpp / lmp;
            }
        }
    }

    for i in 1..RBITS {
        let lmp2 = -(mq.x * mq.x * F12::from_i128(3) / (mq.y * F12::from_i128(2)) * (p.x - mq.x) + mq.y - p.y);
        mq = mq + mq;
        let l2mp =  p.x - mq.x;
        g = g * g * lmp2 / l2mp;
        if bits[i as usize] == 1 {
            if mq.x == q.x  {
                let lmpp = p.x - q.x;                
                mq = q + mq;
                g = g * lmpp ;
            } else {
                let k = (mq.y - q.y) / (mq.x - q.x);
                let lmpp = k * (p.x - q.x) + q.y - p.y;
                mq = q + mq;
                let lmp = p.x - mq.x;
                g = g * lmpp / lmp;
            }
        }
    }
    if RBITSZERO {
        f/g * F12::from_i128(-1)
    } else {
        f/g
    }
}